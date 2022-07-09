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
// for JSON of config data
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/foreach.hpp>
#include <boost/optional.hpp>
using namespace boost::property_tree;


class World {
public:
    string id_;
    const Env env_;
public:
    int priority_mut;
public:
    World(string i, Env e): id_(i), env_(e) {};
    int chronicle(); // main process
    void save_lives();
    void save_labels();
    void save_acts();
};

int World::chronicle() {
    // lives setting
    vector<Lives> lives(env_.C, env_.init_lives);
    LineageLabel labels = env_.init_labels;
    Act act = env_.init_act;
    // recorder
    LivesRec rec_lives(env_.C, env_.L);
    ActRec rec_act(env_.L);
    LineageRec rec_labels(env_.L);
    // recording initial state
    rec_lives.record(lives); rec_act.record(act); rec_labels.record(labels);

    // random generator setting
    std::mt19937 engine = PDS::create_rand_engine();
    std::uniform_real_distribution<double> prob(0, 1);
    std::uniform_real_distribution<double> kh_sample(env_.khs[0], env_.khs[1]);
    std::uniform_real_distribution<double> kp_sample(env_.kps[0], env_.kps[1]);

    // Selection Strength
    int Cs = env_.C * (1.0 - env_.S);

    // Replication & Runge-Kutta system
    runge_kutta4<Lives> RK4;

    // Priority of Mutation
    // Host = 1, Parasite = 2
    priority_mut = 1;

    // Main loop
    for (int rd=0; rd<env_.R; rd++) {
        std::cout << "Round " << rd+1 << std::endl;

        if (prob(engine) < 3/23) {
            priority_mut = 2;
        }
        
        // Selection
        vector<int> killed;
        killed = PDS::selection_mapper(env_.C, Cs);
        for (int k : killed) { for (int l=0; l<env_.L; l++) {
            lives[k][l] = 0.0;
        } }
        std::cout << "Section DONE." << rd+1 << std::endl;

        // Fusion-Division
        vector<vector<int>> fusdived;
        fusdived = PDS::fusdiv_mapper(env_.C, env_.FD);
        for (vector<int> fd : fusdived) {
            vector<int> fused = PDS::sum_2vecs_int(
                lives[fd[0]], lives[fd[1]], env_.L
            );
            // detection empty fusion
            for (int l=0; l<env_.L; l++) {
                if (fused[l] > 0.0) {
                    std::binomial_distribution<> binomial((int)fused[l], 0.5);
                    int oneside = binomial(engine);
                    
                    lives[fd[0]][l] = (double)oneside;
                    lives[fd[1]][l] = (double)(fused[l] - oneside);
                }
            }
        } // Fusion-Division
        std::cout << "Fusion-Division DONE." << rd+1 << std::endl;

        // Decting Extinction
        vector<double> sum_line = PDS::sum_lineage(lives, env_.C, env_.L);
        for (int l=0; l<env_.L; l++) {
            if (sum_line[l] <= 0.0 & labels[l] != 0) {
                // Updating labels
                labels.update_none(l);
                // Updating activity
                PDS::zero_act(act, l, env_.L);
            }
        }
        std::cout << "Extinction Detection DONE." << rd+1 << std::endl;

        //##### Replication #####//
        vector<int> index_host = labels.index_hosts();
        vector<int> index_para = labels.index_parasites();
        vector<int> dev_hosts(env_.C, 0), dev_paras(env_.C, 0);
        Replication rep(act, env_.N, env_.L);
        #pragma omp parallel for
        for (int c=0; c<env_.C; c++) {

        // ### Main Process of Rep ### ///
        double sum_lives = PDS::sum(lives[c], env_.L);
        // detecting host
        int flag_host = 0;
        for (int i : index_host) { if (lives[c][i] > 0) {flag_host += 1;} }
        if (0 < sum_lives & sum_lives < env_.N & flag_host > 0) {
            // Getting initical concs
            for (int i : index_host) { dev_hosts[c] = int(lives[c][i]); }
            for (int i : index_para) { dev_paras[c] = int(lives[c][i]); }
            // Main development
            Lives clives = lives[c];
            boost::numeric::odeint::integrate_const(
                RK4, rep, clives, 0.0, 15.0, 0.01
            );
            lives[c] = clives;
            // Count replication
            for (int i : index_host) {
                dev_hosts[c] = int(lives[c][i]) - dev_hosts[c];
            }
            for (int i : index_para) {
                dev_paras[c] = int(lives[c][i]) - dev_paras[c];
            }
        // ### Main Process of Rep; END ### ///

        } } // Replication
        std::cout << "Replication DONE." << rd+1 << std::endl;

        //##### Mutation #####//
        bool flag_host_mut = false, flag_para_mut = false;
        int index_host_mut = -1, index_para_mut = -1;
        if (labels.detect_none()){
        // ### Mutation by Replication ### ///
        vector<int> index_none = labels.index_nones();
        if (priority_mut == 1) {
        // ## Generation Priority: Parasite < Host ## //
        // Host to Host
        for (int c=0; c<env_.C; c++) {
            if (dev_hosts[c] > 0.0) {
                bool mut_h = (dev_hosts[c] * env_.mhh >= prob(engine));
                if (mut_h) {
                    flag_host_mut = true; lives[c][index_none[0]] += 1.0;
                    index_host_mut = index_none[0];
                    // Updating label
                    labels.update_host(index_host_mut);
        } } }

        // Host to Parasite
        for (int c=0; c<env_.C; c++) {
            if (flag_host_mut && index_none.size() < 2) { break; }
            if (dev_hosts[c] > 0.0) {
                bool mut_p = (dev_hosts[c] * env_.mhp >= prob(engine));
                if (mut_p) {
                    if (flag_para_mut) {
                        lives[c][index_para_mut] += 1.0;
                        // Updating label
                        labels.update_parasite(index_para_mut);
                    } else if (!flag_para_mut && flag_host_mut && index_none.size() > 1) {
                        flag_para_mut = true; lives[c][index_none[1]] += 1.0;
                        index_para_mut = index_none[1];
                        // Updating label
                        labels.update_parasite(index_para_mut);
                    } else if (!flag_host_mut && !flag_para_mut) {
                        flag_para_mut = true; lives[c][index_none[0]] += 1.0;
                        index_para_mut = index_none[0];
                        // Updating label
                        labels.update_parasite(index_para_mut);
                } }
            }
        }
        // Parasite to Parasite
        for (int c=0; c<env_.C; c++) {
            if (flag_host_mut && index_none.size() < 2) { break; }
            if (dev_paras[c] > 0.0) {
                bool mut_p = (dev_paras[c] * env_.mpp >= prob(engine));
                if (mut_p) {
                    if (flag_para_mut) {
                        lives[c][index_para_mut] += 1.0;
                        // Updating label
                        labels.update_parasite(index_para_mut);
                    } else if (!flag_para_mut && flag_host_mut && index_none.size() > 1) {
                        flag_para_mut = true; lives[c][index_none[1]] += 1.0;
                        index_para_mut = index_none[1];
                        // Updating label
                        labels.update_parasite(index_para_mut);
                    } else if (!flag_host_mut && !flag_para_mut) {
                        flag_para_mut = true; lives[c][index_none[0]] += 1.0;
                        index_para_mut = index_none[0];
                        // Updating label
                        labels.update_parasite(index_para_mut);
                } }
            }
        } } else if (priority_mut == 2) {
        // ## Generation Priority: Parasite > Host ## //
        // Host to Parasite
        for (int c=0; c<env_.C; c++) {
            if (dev_hosts[c] > 0.0) {
                bool mut_p = (dev_hosts[c] * env_.mhp >= prob(engine));
                if (mut_p) {
                    if (flag_para_mut) {
                        lives[c][index_para_mut] += 1.0;
                        // Updating label
                        labels.update_parasite(index_para_mut);
                    } else {
                        flag_para_mut = true; lives[c][index_none[0]] += 1.0;
                        index_para_mut = index_none[0];
                        // Updating label
                        labels.update_parasite(index_para_mut);
                    } 
                }
            }
        }
        // Parasite to Parasite
        for (int c=0; c<env_.C; c++) {
            if (dev_paras[c] > 0.0) {
                bool mut_p = (dev_paras[c] * env_.mpp >= prob(engine));
                if (mut_p) {
                    if (flag_para_mut) {
                        lives[c][index_para_mut] += 1.0;
                        // Updating label
                        labels.update_parasite(index_para_mut);
                    } else {
                        flag_para_mut = true; lives[c][index_none[0]] += 1.0;
                        index_para_mut = index_none[0];
                        // Updating label
                        labels.update_parasite(index_para_mut);
                    } 
                }
            }
        }
        // Host to Host
        for (int c=0; c<env_.C; c++) {
            if (flag_para_mut && index_none.size() < 2) { break; }
            if (dev_hosts[c] > 0.0) {
                bool mut_h = (dev_hosts[c] * env_.mhh >= prob(engine));
                if (mut_h) {
                    if (flag_host_mut) {
                        flag_host_mut = true; lives[c][index_none[1]] += 1.0;
                        index_host_mut = index_none[1];
                        // Updating label
                        labels.update_host(index_host_mut);
                    }  else if (!flag_para_mut) {
                        flag_host_mut = true; lives[c][index_none[0]] += 1.0;
                        index_host_mut = index_none[0];
                        priority_mut = 2;
                        // Updating label
                        labels.update_host(index_host_mut);
                    } else if (!flag_host_mut && flag_para_mut && index_none.size() > 1) {
                        flag_host_mut = true; lives[c][index_none[1]] += 1.0;
                        index_host_mut = index_none[1];
                        // Updating label
                        labels.update_host(index_host_mut);
                    }
                }
            }
        }

        } 
        // ### Mutation by Replication; END ### //
        std::cout << "Mutation DONE." << rd+1 << std::endl;

        }

        // ### Updating activity ### //
        // Parasite activity
        if (flag_para_mut > 0) {
            // Parasite replicated
            for (int l=0; l<env_.L; l++) {
                act[index_para_mut][l] = PDS::activity_sampler(
                    kh_sample, kp_sample,
                    labels[index_para_mut], labels[l]
                );
            }
            // Parasite replicating
            for (int l=0; l<env_.L; l++) {
                act[l][index_para_mut] = PDS::activity_sampler(
                    kh_sample, kp_sample,
                    labels[l], labels[index_para_mut]
                );
            }
        }
        // Host activity
        if (flag_host_mut > 0) {
            // Host replicated
            for (int l=0; l<env_.L; l++) {
                act[index_host_mut][l] = PDS::activity_sampler(
                    kh_sample, kp_sample,
                    labels[index_host_mut], labels[l]
                );
            }
            // Host replicating
            for (int l=0; l<env_.L; l++) {
                act[l][index_host_mut] = PDS::activity_sampler(
                    kh_sample, kp_sample,
                    labels[l], labels[index_host_mut]
                );
            }
        }
        // ### Updating activity; END ### //
        std::cout << "Updating Activity DONE." << rd+1 << std::endl;


        // recording conc
        rec_lives.record(lives); rec_act.record(act); rec_labels.record(labels);

        // Display
        PDS::display_lives(lives, env_.C, env_.L);
        PDS::display_label(labels, env_.L);
        PDS::display_act(act, env_.L);

        // If All Hosts extinct, process is broken.
        sum_line = PDS::sum_lineage(lives, env_.C, env_.L);
        int flag_host = 0;
        for (int i : index_host) {
            if (sum_line[i] > 0) {flag_host += 1;}
        }
        if (flag_host == 0) {
            rec_lives.save(id_); rec_act.save(id_); rec_labels.save(id_);
            return 0;
        }
        
    } // Main loop

    rec_lives.save(id_); rec_act.save(id_); rec_labels.save(id_);
    return 1;
} // All hosts extinction: 0, al least one host lineage survival: 1


#endif  // WORLD_HPP_
