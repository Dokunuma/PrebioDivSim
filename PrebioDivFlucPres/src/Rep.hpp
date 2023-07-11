#ifndef REP_HPP_
#define REP_HPP_

#include "Types.hpp"

struct Replication {
public:
    Act act; double N; int L;
public:
    Replication(Act a, double n, int l): act(a), N(n), L(l) {}
    void operator()(const VectorXd& x, VectorXd& dx, double t) {
        // density effect term
        double capa = 1.0 - (x.sum()/N);
        // interaction and development
        VectorXd inter(L); inter.setZero();
        inter = act.values()*x;
        dx = inter * capa;
    }
}; // class Replication

#endif  // REP_HPP_
