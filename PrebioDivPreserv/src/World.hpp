#ifndef WORLD_HPP_
#define WORLD_HPP_

#include "Utils.hpp"
#include "Types.hpp"
#include "Rep.hpp"

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
#include <stdio.h>
#include <fstream>


Lives replication(Lives clives, Act act, double N, int L) {
    Lives vec = clives;
    // Replication & Runge-Kutta system
    runge_kutta4<Lives> RK4;
    Replication rep(act, N, L);

    boost::numeric::odeint::integrate_const(
        RK4, rep, vec, 0.0, 15.0, 0.01
    );

    return vec;
}


class World {
public:
    string id_;
    const Env env_; const Act act_;
public:
    vector<Lives> lives; LivesRec rec;
public:
    World(Env e, Act a): id_("TEST"), env_(e), act_(a) {};
    World(string i, Env e, Act a): id_(i), env_(e), act_(a) {};
    void set_ID(string i) { id_ = i; };
    int chronicle(); // main process
    void save_env();
    void save_act();
};

int World::chronicle() {
    // lives setting
    lives = vector<Lives>(env_.C, env_.init_lives);
    // recorder
    rec = LivesRec(env_.C, env_.L);
    // recording initial conc
    rec.record(lives);

    // the num of host lineages
    int host_num = env_.index_host.size();

    // random generator setting
    std::mt19937 engine = PDS::create_rand_engine();

    // Selection Strength
    int Cs = env_.C * (1.0 - env_.S);

    // Replication & Runge-Kutta system
    Replication rep(act_, env_.N, env_.L);
    runge_kutta4<Lives> RK4;

    // Main loop
    for (int rd=0; rd<env_.R; rd++) {
        std::cout << "Process " << id_ << " Round " << rd+1 << std::endl;
        
        // Selection
        vector<int> killed;
        killed = PDS::selection_mapper(env_.C, Cs, engine);
        for (int k : killed) { for (int l=0; l<env_.L; l++) {
            lives[k][l] = 0.0;
        } }

        // Fusion-Division
        vector<vector<int>> fusdived;
        fusdived = PDS::fusdiv_mapper(env_.C, env_.FD, engine);
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


        // Replication
        #pragma omp parallel for
        for (int c=0; c<env_.C; c++) {
            double sum = PDS::sum(lives[c], env_.L);
            // detecting host
            int flag_host = 0;
            for (int i : env_.index_host) {
                if (lives[c][i] > 0) {flag_host += 1;}
            }
            if (sum > 0 & sum < env_.N & flag_host > 0) {
                Lives clives = replication(lives[c], act_, env_.N, env_.L);
                boost::numeric::odeint::integrate_const(
                    RK4, rep, clives, 0.0, 15.0, 0.01
                );
                lives[c] = clives;
        }   } // Replication

        // recording conc
        rec.record(lives);

        // If Hosts extinct, process is broken.
        Lives sum_line = PDS::sum_lineage(lives, env_.C, env_.L);
        int flag_host = 0;
        for (int i : env_.index_host) {
            if (sum_line[i] > 0) {flag_host += 1;}
        }
        if (flag_host < host_num) {
            std::cout << "Extinction " << id_ << std::endl;
            return 0;
        }
        
    } // Main loop

    std::cout << "Survival " << id_ << std::endl;
    return 1;
}

void World::save_env() {
    /// save path
    std::string results, ext;
    results = string("./results/env/"); ext = string(".json");
    string path = results + id_ + ext;
    
    /// constructing config
    ptree config;
    config.put("Env.Round", env_.R);
    config.put("Env.Compartment", env_.C);
    config.put("Env.Selection", env_.S);
    config.put("Env.Fusion_Division", env_.FD);
    config.put("Env.Capacity", env_.N);
    config.put("Env.Lineage", env_.L);
    // initial conc
    ptree iconc;
    for (int l=0; l<env_.L; l++) {
        string index = std::to_string(l+1);
        iconc.put(index, env_.init_lives[l]);
    }
    config.add_child("Env.Init_Lives", iconc);

    /// Writing
    write_json(path, config);
}

void World::save_act() {
    std::cout << "Activity saver " << id_ << std::endl;
    // Saver
    std::ofstream writing_file;
    std::string filename = "./results/act/" + id_ + ".csv";
    writing_file.open(filename, std::ios::out);

    // Main
    for (int l1=0; l1<env_.L; l1++) { 
    for (int l2=0; l2<env_.L; l2++) {
        if (l1 == env_.L-1) {
            if (l2 == env_.L-1) {
                writing_file << act_[l1][l2] << std::endl;
            } else {
                writing_file << act_[l1][l2] << ",";
            }
        } else {
            writing_file << act_[l1][l2] << ",";
        }
    } }
    writing_file.close();
    std::cout << "Activity save DONE " << id_ << std::endl;
}


#endif  // WORLD_HPP_
