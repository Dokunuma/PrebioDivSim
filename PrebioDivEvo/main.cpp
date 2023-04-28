#include "src/World.hpp"
#include "src/Utils.hpp"
#include "src/Types.hpp"
#include <sys/stat.h>
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
const array<double, 2> range_kh = {1.0, 3.0};
const array<double, 2> range_kp = {0.0, 10.0};
// Mutation rates
const double mhh = 0.02, mpp = 0.002, mhp = 0.001;
// Trial number in a process
const int TRIAL = 2;


int create(string code) {
    std::cout << "Process " << code << std::endl;
    
    vector<double> init_lives(L);
    init_lives[0] = N;

    vector<int> init_lineages(L, 0);
    init_lineages[0] = 1;

    MatrixXd init_acts = MatrixXd::Constant(L, L, 0.0);
    init_acts(0, 0) = 2.0;

    Env env(
        R, C, S, FD, N, L, init_lives, init_acts, init_lineages, range_kh, range_kp,
        mhh, mpp, mhp, code
    );
    World world(code, env);

    int div = world.chronicle();

    return div;
}


int main() {
    for (int i=0; i<TRIAL; i++) {
        string code = PDS::random_code_generator(10);
        //　getting time and covert it to string
        string now = PDS::get_time_str();
        string code_exp = now + "_" + code;
        create(code_exp);
    }

    return 0;
}
