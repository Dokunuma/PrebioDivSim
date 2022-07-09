#ifndef REP_HPP_
#define REP_HPP_

#include "Types.hpp"

struct Replication {
public:
    Act act; double N; int L;
public:
    Replication(Act a, double n, int l): act(a), N(n), L(l) {}
    void operator()(const Lives& x, Lives& dx, double t) {
        // density effect term
        double capa = 1.0 - (PDS::sum(x, L)/N);
        // interaction and development
        for (int i=0; i<L; i++) { double inter = 0.0;
            for (int j=0; j<L; j++) {
                inter += act[i][j] * x[j];
            }
            dx[i] = x[i] * inter * capa;
    } }
}; // class Replication

#endif  // REP_HPP_
