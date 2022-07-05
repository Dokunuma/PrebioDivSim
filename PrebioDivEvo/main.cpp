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

#include <fstream>

#include <cassert>

// Environmental arameters
const int R = 1000, C = 3000, FD = 5000, L = 5;
const double S = 0.25, N = 100.0;
// Ranges of activities of host and parasite replications
const vector<double> range_kh = {1.0, 3.0};
const vector<double> range_kp = {0.0, 10.0};
// Mutation rates
const double mhh = 0.02, mpp = 0.002, mhp = 0.001;
// Trial number in a process
const int TRIAL = 2;


int create(Act act, string code) {
    std::cout << "Process " << code << std::endl;
    Lives ilives(L, 0.0); ilives[0] = N;
    LineageLabel ilabels(vector<int>(L, 0));
    ilabels.update_host(0);
    Env env(
        R, C, S, FD, N, ilives, L, range_kh, range_kp, ilabels, act,
        mhh, mpp, mhp
    );
    World world(code, env);

    int div = world.chronicle();

    return div;
}


int main() {
    string code = PDS::random_code_generator(10);
    std::cout << "Process Start: " << code << std::endl;

    vector<int> div = {};

    int n = 0;
    Act act = PDS::zero_act_generator(L);
    act[0][0] = 2.0;

    for (int i=0; i<TRIAL; i++) {
        create(act, code+"_"+std::to_string(i));
    }

    return 0;
}
