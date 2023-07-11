#ifndef TYPES_HPP_
#define TYPES_HPP_

#if __INTELLISENSE__
#undef __ARM_NEON
#undef __ARM_NEON__
#endif

#include <iostream>
#include <fstream>
#include <array>
#include <vector>
#include <string>
#include <cassert>
#include <Eigen/Dense>
#include <nlohmann/json.hpp>

using std::ofstream;
using std::array;
using std::vector;
using std::string;

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::VectorXi;


class Live {
    vector<VectorXd> values_;
    int comp_size_;
    int lineage_size_;

public:
    Live(int, int);
    double value(int i, int j) { return values_[i](j); }
    VectorXd values(int i) { return values_[i]; };
    void filled(double);
    void filled(vector<double>);
    void update(int i, int j, double v) { values_[i](j) = v; }
    void update(int i, double v);
    void update(int i, VectorXd v) { values_[i] = v; }
    void save(ofstream *ofs);
    VectorXd sum_lineage();
    VectorXd sum_compartment();
    double sum(int, vector<int>);
    double sum_compartment(int i) { return values_[i].sum(); };
    bool exist(int, vector<int>);
    void display_lineage();
    VectorXi fuse(int, int);
};

Live::Live(int c, int l) {
    comp_size_ = c; lineage_size_ = l;
    VectorXd vec(lineage_size_); vec.setZero();
    for (int i=0; i<c; i++) {
        values_.push_back(vec);
    }
}

void Live::filled(double v) {
    for (int i=0; i<comp_size_; i++) { for (int j=0; j<lineage_size_; j++) {
        values_[i](j) = v;
    }}
}

void Live::filled(vector<double> vec) {
    for (int i=0; i<comp_size_; i++) { for (int j=0; j<lineage_size_; j++) {
        values_[i](j) = vec[j];
    }}
}

inline void Live::update(int c, double v) {
    for (int l=0; l<lineage_size_; l++) { values_[c](l) = v; }
}

void Live::save(ofstream *ofs) {
    for (int i=0; i<comp_size_; i++) { for (int j=0; j<lineage_size_; j++) {
        if (i < comp_size_-1 || j < lineage_size_-1) {
            *ofs << values_[i](j) << ",";
        } else {
            *ofs << values_[i](j) << std::endl;
        }
    }}
}

VectorXd Live::sum_lineage() {
    VectorXd s(lineage_size_); s.setZero();

    for (int i=0; i<comp_size_; i++) { for (int j=0; j<lineage_size_; j++) {
        s(j) += values_[i](j);
    }}

    return s;
}

VectorXd Live::sum_compartment() {
    VectorXd s(comp_size_); s.setZero();

    for (int i=0; i<comp_size_; i++) { for (int j=0; j<lineage_size_; j++) {
        s(i) += values_[i](j);
    }}

    return s;
}

inline double Live::sum(int c, vector<int> index) {
    double ans = 0.0;
    for (int i : index) { ans += values_[c](i); }
    return ans;
}

inline bool Live::exist(int c, vector<int> index) {
    for (int i : index) { if (values_[c](i) > 0.0) { return true; } }
    return false;
}

void Live::display_lineage() {
    VectorXd l = sum_lineage();
    for (int i=0; i<lineage_size_; i++){
        if (i < lineage_size_-1) { std::cout << l[i] << ", "; }
        else { std::cout << l[i] << std::endl; }
    }
}

inline VectorXi Live::fuse(int i, int j) {
    VectorXi vec = (values_[i]+values_[j]).cast <int> ();
    return vec;
};


class Act {
    MatrixXd values_;
    int size_;

public:
    Act(int n) {values_.resize(n, n); size_=n; initilized();}
    void initilized();
    MatrixXd values() {return values_;}
    void update(int i, int j, double v) {values_(i, j) = v;}
    void update(MatrixXd);
    void save(ofstream *ofs);
    void zero(int);
};

void Act::initilized() {
    for (int i=0; i<size_; i++) { for (int j=0; j<size_; j++) {
        values_(i, j) = 0.0;
    }}
}

void Act::update(MatrixXd m) {
    // Comparing matirices's size
    assert(size_ == m.rows()); assert(size_ == m.cols());
    values_ = m;
}

void Act::save(ofstream *ofs) {
    for (int i=0; i<size_; i++) { for (int j=0; j<size_; j++) {
        if (i < size_-1 || j < size_-1) {
            *ofs << values_(i, j) << ",";
        } else {
            *ofs << values_(i, j) << std::endl;
        }
    }}
}

void Act::zero(int l) {
    for (int i=0; i<size_; i++) {
        values_(l, i) = 0.0; values_(i, l) = 0.0;
    }
}


class LineageLabel {
    // None = 0, Host = 1, Parasite = 2
    vector<int> values_;
    int size_;

public:
    LineageLabel(int s) {values_.reserve(s); size_ = s;};
    int operator[](int i) {return values_[i];}
    int size() {return size_;}
    void update(vector<int> v) {values_ = v;}
    void update_none(int i) {values_[i] = 0;}
    void update_host(int i) {values_[i] = 1;}
    void update_parasite(int i) {values_[i] = 2;}
    void update_para(int i) {update_parasite(i);}
    bool detect_none();
    bool detect_host();
    bool detect_parasite();
    bool detect_para() { return detect_parasite(); }
    vector<int> index_nones();
    vector<int> index_hosts();
    vector<int> index_parasites();
    vector<int> index_para() {return index_parasites();}
    void save(ofstream *ofs);
};

inline bool LineageLabel::detect_none() {
    for (int v : values_) { if (v == 0) {return true;} }
    return false;
}

inline bool LineageLabel::detect_host() {
    for (int v : values_) { if (v == 1) {return true;} }
    return false;
}

inline bool LineageLabel::detect_parasite() {
    for (int v : values_) { if (v == 2) {return true;} }
    return false;
}

inline vector<int> LineageLabel::index_nones() {
    vector<int> index = {};
    for (int i=0; i<size_; i++) {
        if (values_[i] == 0) {index.push_back(i);}
    }
    return index;
}

inline vector<int> LineageLabel::index_hosts() {
    vector<int> index = {};
    for (int i=0; i<size_; i++) {
        if (values_[i] == 1) {index.push_back(i);}
    }
    return index;
}

inline vector<int> LineageLabel::index_parasites() {
    vector<int> index = {};
    for (int i=0; i<size_; i++) {
        if (values_[i] == 2) {index.push_back(i);}
    }
    return index;
}

inline void LineageLabel::save(ofstream *ofs) {
    for (int i=0; i<size_; i++) {
        if (i < size_-1) { *ofs << values_[i] << ","; }
        else { *ofs << values_[i] << std::endl; }
    }
}


class Env {
public:
    int R_, C_, FD_, L_; double S_, N_; 
    vector<double> init_live_;
    MatrixXd init_act_; vector<int> init_lineage_;
    array<double, 2> kh_range_, kp_range_;
    double mhh_, mpp_, mhp_;
    string id_;

public:
    Env(
        int r, int c, double s, int fd, double n, int l,
        vector<double> init_live, MatrixXd init_act, vector<int> init_lineage,
        array<double, 2> kh, array<double, 2> kp,
        double mhh, double mpp, double mhp, string id
    ):
        R_(r), C_(c), S_(s), FD_(fd), N_(n), L_(l), 
        init_live_(init_live), init_act_(init_act), init_lineage_(init_lineage),
        kh_range_(kh), kp_range_(kp), mhh_(mhh_), mpp_(mpp), mhp_(mhp), id_(id) {};
    void save();
};

void Env::save() {
    nlohmann::json j;

    j["round"] = R_; j["compartment"] = C_; j["fusdiv"] = FD_; j["lineage"] = L_;
    j["selection"] = S_; j["nutrient"] = N_;
    j["kh"] = kh_range_;j["kp"] = kp_range_;
    j["h_to_h"] = mhh_; j["p_to_p"] = mpp_; j["h_to_p"] = mhp_;
    j["ID"] = id_;

    string dir;
    dir = "./results/envs/" + id_ + ".json";

    ofstream f(dir);
    f << std::setw(4) << j << std::endl;
}

#endif  // TYPES_HPP_
