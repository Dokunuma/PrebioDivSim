#ifndef TYPES_HPP_
#define TYPES_HPP_

#include <iostream>
#include <fstream>
#include <array>
#include <vector>
#include <string>
#include <cassert>

using std::array;
using std::vector;
using std::string;

// live conc in each compartment
using Lives = vector<double>;
// activity of replication
using Act = vector<vector<double>>;


namespace PDS {
Lives sum_lineage(vector<Lives> lives, int compartment, int lineage);
}


class Env {
public:
    const int R, C, FD, L; const double S, N; 
    const Lives init_lives;
    const vector<int> index_host;
public:
    Env(): R(0), C(0), S(0.0), FD(0), N(0.0), init_lives(Lives{0.0}), L(1) {};
    Env(
        int r, int c, double s, int fd, double n, Lives il, int l,
        vector<int> ih
    ):
        R(r), C(c), S(s), FD(fd), N(n), init_lives(il), L(l), index_host(ih) {
            assert(L == init_lives.size());
    }; // Env's constructor
}; // class Env


class LivesRec {
public:
    vector<Lives> data;
    int size, compartment, lineage;

    LivesRec() {};
    LivesRec(int c, int l): compartment(c), lineage(l) { data = {}; size = 0; }
    void record(vector<Lives>);
    void save(string);
    void display();
};

inline void LivesRec::record(vector<Lives> lives) {
    Lives sum = PDS::sum_lineage(lives, compartment, lineage);

    size += 1;
    data.push_back(sum);
}

inline void LivesRec::save(string id) {
    std::string results, ext;
    results = string("./results/pop/"); ext = string(".csv");
    string path = results + id + ext;
    std::ofstream ofs(path);

    for (int s=0; s<size; s++) { for (int l=0; l<lineage; l++) {
        if (l+1!=lineage) { ofs << data[s][l] << ","; }
        else { ofs << data[s][l] << std::endl; }
    } }
}

inline void LivesRec::display() {
    for (int i=0; i<size; i++) { 
        std::cout << "{ ";
        for (int l=0; l<lineage; l++) { std::cout << data[i][l] << " "; }
        std::cout << "}" << std::endl;
    }
}


namespace PDS {

inline double sum(Lives lives, int lineage) {
    double sum = 0.0;
    for (int l=0; l<lineage; l++) {
        sum += lives[l];
    }

    return sum;
}

vector<double> sum_2Lives_d(Lives v1, Lives v2, int dim) {
    vector<double> sum(dim, 0.0);
    for (int i=0; i<dim; i++) { sum[i] = v1[i] + v2[i]; }

    return sum;
}

vector<int> sum_2Lives_i(Lives v1, Lives v2, int dim) {
    vector<int> sum(dim, 0.0);
    for (int i=0; i<dim; i++) { sum[i] = (int)(v1[i] + v2[i]); }

    return sum;
}

vector<int> sum_2vecs_int(vector<double> v1, vector<double> v2, int dim) {
    vector<int> sum(dim, 0.0);
    for (int i=0; i<dim; i++) { sum[i] = (int)(v1[i] + v2[i]); }

    return sum;
}

Lives sum_lineage(vector<Lives> lives, int compartment, int lineage) {
    Lives sum(lineage, 0.0);

    for (int c=0; c<compartment; c++) { for (int l=0; l<lineage; l++) {
        sum[l] += lives[c][l];
    } }

    return sum;
}

Lives sum_compartment(vector<Lives> lives, int compartment, int lineage) {
    Lives sum(compartment, 0.0);

    for (int c=0; c<compartment; c++) { for (int l=0; l<lineage; l++) {
        sum[c] += lives[c][l];
    } }

    return sum;
}

}

#endif  // TYPES_HPP_
