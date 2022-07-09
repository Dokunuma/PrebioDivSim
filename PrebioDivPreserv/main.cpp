#include "src/World.hpp"
#include "src/Utils.hpp"
#include "src/Types.hpp"
#include <sys/stat.h>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/foreach.hpp>
#include <boost/optional.hpp>
#include <random>

#include <iostream>
#include <iterator>
#include <algorithm>

#include <stdio.h>
#include <fstream>
#include <string>

#include <cassert>

// Environmental parameters
const int R = 100, C = 3000, FD = 5000, L = 2;
const double S = 0.25, N = 100.0;
// Index for host lineage
const vector<int> index_host = {0};
// Replication activity values of host and parasite
const vector<double> khs = {0.2, 1.7, 2.0, 2.3, 2.6};
const vector<double> kps = {0.2, 1.7, 2.0, 2.3, 2.6, 7.0, 10, 20};
// Trial number of each condition
const int TRIAL = 100;


int creat_each(Act act, int id) {
    string code = PDS::random_code_generator(10) + "_" + std::to_string(id);

    std::cout << "Process " << code << std::endl;
    Lives ilives(2, N/2);
    Env env(R, C, S, FD, N, ilives, L, index_host);
    World world(code, env, act);
    int div = world.chronicle();
    world.save_act();
    std::cout << "Process DONE " << code << std::endl;

    return div;
}

vector<int> create(Act act, string code) {
    vector<int> divs = {};

    for (int n=0; n<TRIAL; n++) {
        std::cout << "Process " << code << " #" << n << std::endl;
        Lives ilives(L, N/L);
        Env env(R, C, S, FD, N, ilives, L, index_host);
        World world(code+"-"+std::to_string(n), env, act);
        int div = world.chronicle();
        divs.push_back(div);

        if (n==TRIAL-1) { world.save_act(); }
    }
    return divs;
}


vector<Act> create_acts(int start, int end) {
    vector<Act> acts = {};
    int i = 0;
    for (double kh : khs) {
    for (double kp : kps) {
        if (start <= i & i <= end & kh <= kp) {
            Act act{{
                {kh, 0.0}, {kp, 0.0}
            }};
            acts.push_back(act);
        }
        i += 1;
    } }

    return acts;
}

void vector_saver(vector<int> vec, int size, string code) {
    std::cout << "Result Saver " << code << std::endl;
    // Saver
    std::ofstream writing_file;
    std::string filename = "./results/pop/" + code + ".csv";
    writing_file.open(filename, std::ios::out);

    // Main
    for (int i=0; i<size; i++) {
        if (i == size-1) {
            writing_file << vec[i] << std::endl;
        } else {
            writing_file << vec[i] << ",";
        }
    }
    writing_file.close();
    std::cout << "Result Save DONE " << code << std::endl;
}

int main(int argc, char *argv[]) {
    // Experimental Code
    std::cout << "##### Process Start #####" << std::endl << std::endl;

    string START_str = argv[1], END_str = argv[2];
    int START = std::atoi(START_str.c_str()), END = std::atoi(END_str.c_str());

    vector<Act> acts = create_acts(START, END);
    int n = 0;
    for (Act act : acts) {
        string code = PDS::random_code_generator(10);
        vector<int> divs  = create(act, code);

        std::cout << "Result { ";
        for (int i=0; i<TRIAL; i++) { std::cout << divs[i] << " "; }
        std::cout << "}" << std::endl;

        vector_saver(divs, TRIAL, code + "-" + std::to_string(n));
        std::cout << std::endl;
        n += 1;
    }

    return 0;
}
