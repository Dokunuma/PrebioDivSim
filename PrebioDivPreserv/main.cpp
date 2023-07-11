#include "src/World.hpp"
#include "src/Utils.hpp"
#include "src/Types.hpp"
#include <sys/stat.h>
#include <random>

#include <iostream>
#include <iterator>
#include <algorithm>

#include <stdio.h>
#include <fstream>
#include <string>

#include <cassert>

// Environmental parameters
const int R = 100, C = 3000, FD = 5000, L = 3;
const double S = 0.25, N = 100.0;
// Index for host lineage
const vector<int> index_host = {0};
// Trial number of each condition
const int TRIAL = 100;


int create() {
    string code = PDS::random_code_generator(10);

    MatrixXd act(L, L);
    act << 2.0, 2.0, 0.0, 2.0, 2.0, 0.0, 7.0, 0.1, 0.0;
    Act acts(L); acts.update(act);

    std::cout << "Process " << code << std::endl;
    vector<double> init_lives{N, 0.0, 0.0};

    Env env(R, C, S, FD, N, L, init_lives, index_host, code);
    World world(code, env, acts);
    world.chronicle();

    return 0;
}

int main(int argc, char *argv[]) {
    create();

    return 0;
}
