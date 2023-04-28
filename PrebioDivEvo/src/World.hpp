#ifndef WORLD_HPP_
#define WORLD_HPP_

#include "Utils.hpp"
#include "Types.hpp"
#include "Rep.hpp"
#include "Sample.hpp"

#include <random>
#include <boost/numeric/odeint.hpp>
#include <omp.h>
using namespace boost::numeric::odeint;


class World {
public:
    string id_;
    Env env_;
public:
    World(string i, Env e): id_(i), env_(e) {};
    int chronicle();
};

int World::chronicle() {
    env_.save();
    // lives setting
    Live lives(env_.C_, env_.L_);
    lives.filled(env_.init_live_);
    std::cout << lives.values(0) << std::endl;
    // act matrix setting
    Act acts(env_.L_);
    acts.update(env_.init_act_);
    // lineage label setting
    LineageLabel lineages(env_.L_);
    lineages.update(env_.init_lineage_);

    // output filestream
    ofstream live_ofs("./results/pops/"+env_.id_+".csv");
    ofstream acts_ofs("./results/acts/"+env_.id_+".csv");
    ofstream lins_ofs("./results/lins/"+env_.id_+".csv");

    // saving initial condition
    lives.save(&live_ofs); acts.save(&acts_ofs); lives.save(&live_ofs);

    // random generator setting
    std::mt19937 engine = PDS::create_rand_engine();
    std::uniform_real_distribution<double> prob(0, 1);
    std::uniform_real_distribution<double> kh_sample(env_.kh_range_[0], env_.kh_range_[1]);
    std::uniform_real_distribution<double> kp_sample(env_.kp_range_[0], env_.kp_range_[1]);

    // Selection Strength
    int Cs = env_.C_ * (1.0 - env_.S_);

    // Replication & Runge-Kutta system
    runge_kutta4<VectorXd> RK4;

    // Priority of Mutation
    // Host = 1, Parasite = 2
    int priority_mut = 1;

    // Main loop
    for (int rd=0; rd<env_.R_; rd++) {
        std::cout << "Round " << rd+1 << std::endl;
        
        // Selection
        vector<int> killed = PDS::selection_mapper(env_.C_, Cs);
        for (int k : killed) { lives.update(k, 0.0); }
        std::cout << "Selection DONE. " << rd+1 << std::endl;

        // Fusion-Division
        vector<vector<int>> fusdived;
        fusdived = PDS::fusdiv_mapper(env_.C_, env_.FD_);
        for (vector<int> fd : fusdived) {
            VectorXd fused(env_.L_); fused = lives.fuse(fd[0], fd[1]);
            std::cout << "Fusion-Division DONE." << rd+1 << std::endl;
            // detection empty fusion
            for (int l=0; l<env_.L_; l++) {
                if (fused[l] > 0.0) {
                    std::binomial_distribution<> binomial((int)fused[l], 0.5);
                    double oneside = binomial(engine);
                    
                    lives.update(fd[0], l, oneside);
                    lives.update(fd[1], l, fused[l]-oneside);
                }
            }
        } // Fusion-Division
        std::cout << "Fusion-Division DONE." << rd+1 << std::endl;

        // Decting Extinction
        VectorXd sum_line = lives.sum_lineage();
        for (int l=0; l<env_.L_; l++) {
            if (sum_line[l] <= 0.0 & lineages[l] != 0) {
                // Updating labels
                lineages.update_none(l);
                // Updating activity
                acts.zero(l);
            }
        }
        std::cout << "Extinction Detection DONE." << rd+1 << std::endl;

        //##### Replication #####//
        vector<int> index_host = lineages.index_hosts();
        vector<int> index_para = lineages.index_para();
        vector<int> dev_hosts(env_.C_, 0), dev_paras(env_.C_, 0);
        Replication rep(acts, env_.N_, env_.L_);
        #pragma omp parallel for
        for (int c=0; c<env_.C_; c++) {

        // ### Main Process of Rep ### ///
        double sum_lives = lives.sum_compartment(c);
        // detecting host
        int flag_host = 0;
        if (0 < sum_lives & sum_lives < env_.N_ & lives.exist(c, index_host)) {
            // Getting initical concs
            dev_hosts[c] = int(lives.sum(c, index_host));
            dev_paras[c] = int(lives.sum(c, index_para));
            // Main development
            VectorXd clives = lives.values(c);
            boost::numeric::odeint::integrate_const(
                RK4, rep, clives, 0.0, 15.0, 0.01
            );
            lives.update(c, clives);
            // Count replication
            dev_hosts[c] = int(lives.sum(c, index_host)) - dev_hosts[c];
            dev_paras[c] = int(lives.sum(c, index_para)) - dev_paras[c];
        // ### Main Process of Rep; END ### ///

        } } // Replication
        std::cout << "Replication DONE." << rd+1 << std::endl;

        //##### Mutation #####//
        // Mutation priority
        if (prob(engine) < env_.mhh_/(env_.mhh_+env_.mhp_+env_.mpp_))
            { priority_mut = 1; } else { priority_mut = 2; }
        
        vector<int> compartment_dev_host = PDS::where_more(dev_hosts, 0.0);
        vector<int> compartment_dev_para = PDS::where_more(dev_paras, 0.0);
        bool flag_host_mut = false, flag_para_mut = false;
        int index_host_mut = -1, index_para_mut = -1;
        if (lineages.detect_none()){

        // ### Mutation by Replication ### ///
        vector<int> index_none = lineages.index_nones();
        if (priority_mut == 1) {
        // ## Generation Priority: Parasite < Host ## //
        // Host to Host
        for (int c : compartment_dev_host) {
            if (dev_hosts[c] * env_.mhh_ >= prob(engine)) {
                lives.update(c, index_none[0], 1.0);
                if (!flag_host_mut){
                    flag_host_mut = true;
                    index_host_mut = index_none[0];
                    // Updating label
                    lineages.update_host(index_host_mut);
                }
        } }

        // If host mutated and lineage are vacant, para can mutate.
        if (flag_host_mut && index_none.size() >= 2) { index_para_mut = index_none[1]; }
        else if (flag_host_mut && index_none.size() < 2) { continue; }
        else if (!flag_host_mut) { index_para_mut = index_none[0]; }
        else { continue; }

        // Host to Parasite
        for (int c : compartment_dev_host) {

            if (dev_hosts[c] * env_.mhp_ >= prob(engine)) {
                lives.update(c, index_para_mut, 1.0);
                if (!flag_para_mut) {
                    flag_para_mut = true;
                    lineages.update_parasite(index_para_mut);
                }
            }
        }
        // Parasite to Parasite
        for (int c : compartment_dev_para) {
            if (dev_paras[c] * env_.mpp_ >= prob(engine)) {
                lives.update(c, index_para_mut, 1.0);
                if (!flag_para_mut) {
                    flag_para_mut = true;
                    lineages.update_parasite(index_para_mut);
                } 
            }
        } } else if (priority_mut == 2) {
        // ## Generation Priority: Parasite > Host ## //
        // Host to Parasite
        for (int c : compartment_dev_host) {
            if (dev_hosts[c] * env_.mhp_ >= prob(engine)) {
                lives.update(c, index_none[0], 1.0);
                if (!flag_para_mut) {
                    flag_para_mut = true;
                    index_host_mut = index_none[0];
                    // Updating label
                    lineages.update_host(index_host_mut);
                }
            }
        }
        // Parasite to Parasite
        for (int c : compartment_dev_para) {
            if (dev_paras[c] * env_.mpp_ >= prob(engine)) {
                lives.update(c, index_none[0], 1.0);
                if (!flag_para_mut) {
                    flag_para_mut = true;
                    index_host_mut = index_none[0];
                    // Updating label
                    lineages.update_host(index_host_mut);
                }
            }
        }

        // If host mutated and lineage are vacant, para can mutate.
        if (flag_para_mut && index_none.size() >= 2) { index_host_mut = index_none[1]; }
        else if (flag_para_mut && index_none.size() < 2) { continue; }
        else if (!flag_para_mut) { index_host_mut = index_none[0]; }
        else { continue; }

        // Host to Host
        for (int c : compartment_dev_host) {
            if (dev_hosts[c] * env_.mhh_ >= prob(engine)) {
                lives.update(c, index_host_mut, 1.0);
                if (!flag_host_mut) {
                    flag_host_mut = true;
                    // Updating label
                    lineages.update_host(index_host_mut);
                }
            }
        }  

        } 
        // ### Mutation by Replication; END ### //
        std::cout << "Mutation DONE." << rd+1 << std::endl;

        }

        // ### Updating activity ### //
        // Parasite activity    
        if (flag_para_mut) {
            for (int l=0; l<env_.L_; l++) {
                // Parasite replicated
                double d1 = PDS::activity_sampler(
                    kh_sample, kp_sample,
                    lineages[index_para_mut], lineages[l]
                );
                acts.update(index_para_mut, l, d1);

                // Parasite replicating
                double d2 = PDS::activity_sampler(
                    kh_sample, kp_sample,
                    lineages[l], lineages[index_para_mut]
                );
                acts.update(l, index_para_mut, d2);
            }
        }
        // Host activity
        if (flag_host_mut) {
            for (int l=0; l<env_.L_; l++) {
                // Host replicated
                double d1 = PDS::activity_sampler(
                    kh_sample, kp_sample,
                    lineages[index_host_mut], lineages[l]
                );
                acts.update(index_host_mut, l, d1);
                // Host replicating
                double d2 = PDS::activity_sampler(
                    kh_sample, kp_sample,
                    lineages[l], lineages[index_host_mut]
                );
                acts.update(l, index_host_mut, d2);
            }
        }
        // ### Updating activity; END ### //
        std::cout << "Updating Activity DONE." << rd+1 << std::endl;

        // saving initial condition
        lives.save(&live_ofs); acts.save(&acts_ofs); lives.save(&live_ofs);

        // Display

        // If All Hosts extinct, process is broken.
        if (!lineages.detect_host()) {
            lives.save(&live_ofs); acts.save(&acts_ofs); lives.save(&live_ofs);
            return 0;        
        }
    }
    lives.save(&live_ofs); acts.save(&acts_ofs); lives.save(&live_ofs);
    return 1;
} // All hosts extinction: 0, al least one host lineage survival: 1


#endif  // WORLD_HPP_
