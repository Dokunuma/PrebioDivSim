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
    Act acts_;
public:
    World(string i, Env e, Act a): id_(i), env_(e), acts_(a) {};
    int chronicle();
};

int World::chronicle() {
    // lives setting
    Live lives(env_.C_, env_.L_);
    lives.filled(env_.init_live_);

    // output filestream
    ofstream live_ofs("./results/pops/"+env_.id_+".csv");
    ofstream acts_ofs("./results/acts/"+env_.id_+".csv");

    // saving initial condition
    lives.save(&live_ofs); acts_.save(&acts_ofs);

    // random generator setting
    std::mt19937 engine = PDS::create_rand_engine();

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
            VectorXi fused(env_.L_); fused = lives.fuse(fd[0], fd[1]);
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

        /// ### Replication ### ///
        Replication rep(acts_, env_.N_, env_.L_);
        //#pragma omp parallel for
        for (int c=0; c<env_.C_; c++) {

        double sum_lives = lives.sum_compartment(c);
        double sum_hosts = lives.sum(c, env_.host_index_);
        // detecting host
        int flag_host = 0;
        if (0.0 < sum_lives & sum_lives < env_.N_ & sum_lives > 0.0) {
            // Main development
            VectorXd clives(env_.L_); clives = lives.values(c);

            boost::numeric::odeint::integrate_const(
                RK4, rep, clives, 0.0, 15.0, 0.01
            );
            lives.update(c, clives);
        // ### Main Process of Rep; END ### ///

        } } // Replication
        std::cout << "Replication DONE." << rd+1 << std::endl;

        // saving initial condition
        lives.save(&live_ofs); acts_.save(&acts_ofs);

        // If All Hosts extinct, process is broken.
        if (!lives.exist(env_.host_index_)) {
            lives.save(&live_ofs); acts_.save(&acts_ofs);
            return 0;        
        }
    }
    lives.save(&live_ofs); acts_.save(&acts_ofs);
    return 1;
} // All hosts extinction: 0, al least one host lineage survival: 1


#endif  // WORLD_HPP_
