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
vector<double> sum_lineage(vector<Lives> lives, int compartment, int lineage);
}


class LivesRec {
public:
    vector<vector<double>> data_;
    int size_, C, L;

    LivesRec() {};
    LivesRec(int c, int l): C(c), L(l) { data_ = {}; size_ = 0; }
    int size() {return size_;}
    void record(vector<Lives>);
    void save(string);
};

inline void LivesRec::record(vector<Lives> lives) {
    vector<double> sum = PDS::sum_lineage(lives, C, L);

    data_.push_back(sum); size_ += 1;
}

inline void LivesRec::save(string id) {
    // Saver
    std::ofstream writing_file;
    std::string filename = "./results/pop/" + id + ".csv";
    writing_file.open(filename, std::ios::out);

    // Main
    for (int r=0; r<size_; r++) { for (int l=0; l<L; l++) {
        if (l == L-1) {
            writing_file << data_[r][l] << std::endl;
        } else {
            writing_file << data_[r][l] << ",";
        }
    } }
    writing_file.close();
}


class ActRec {
public:
    vector<vector<double>> data_;
    const int L;
public:
    ActRec(): L(0) {};
    ActRec(int l): L(l) {data_ = {};};
    void record(Act);
    int size() {return data_.size();}
    int lineage() {return L;}
    void save(string);
};

void ActRec::record(Act act) {
    // Error Detection
    assert(act.size() == L);
    
    vector<double> a = {};
    for (int l1=0; l1<L; l1++) {
        for (int l2=0; l2<L; l2++) {
            a.push_back(act[l1][l2]);
        }
    }
    
    data_.push_back(a);
}

inline void ActRec::save(string id) {
    // Saver
    std::ofstream writing_file;
    std::string filename = "./results/act/" + id + ".csv";
    writing_file.open(filename, std::ios::out);

    // Main
    for (vector<double> a : data_) {
        for (int l=0; l<L*L; l++) {
            if (l == L*L-1) {
                writing_file << a[l] << std::endl;
            } else {
                writing_file << a[l] << ",";
            }
        }
    }
    writing_file.close();
}


class LineageLabel {
    // None = 0, Host = 1, Parasite = 2
public:
    vector<int> labels_;
    int dim_;
public:
    LineageLabel();
    LineageLabel(vector<int> l) {labels_ = l; dim_ = labels_.size();}
    int operator[](int i) {return labels_[i];}
    int dim() {return labels_.size();}
    vector<int> values() {return labels_;}
    void update_none(int i) {labels_[i] = 0;}
    void update_host(int i) {labels_[i] = 1;}
    void update_parasite(int i) {labels_[i] = 2;}
    bool detect_none();
    int index_none();
    vector<int> index_nones();
    vector<int> index_hosts();
    vector<int> index_parasites();
};

inline bool LineageLabel::detect_none() {
    for (int l : labels_) {
        if (l == 0) {return true;}
    }
    return false;
}

inline int LineageLabel::index_none() {
    for (int i=0; i<dim_; i++) {
        if (labels_[i] == 0) {return i;}
    }
    return -1;
}

inline vector<int> LineageLabel::index_nones() {
    vector<int> index = {};
    for (int i=0; i<dim_; i++) {
        if (labels_[i] == 0) {index.push_back(i);}
    }
    return index;
}

inline vector<int> LineageLabel::index_hosts() {
    vector<int> index = {};
    for (int i=0; i<dim_; i++) {
        if (labels_[i] == 1) {index.push_back(i);}
    }
    return index;
}

inline vector<int> LineageLabel::index_parasites() {
    vector<int> index = {};
    for (int i=0; i<dim_; i++) {
        if (labels_[i] == 2) {index.push_back(i);}
    }
    return index;
}


class LineageRec {
public:
    vector<vector<int>> lineages_;
    const int dim_;
public:
    LineageRec(): dim_(0) {};
    LineageRec(int d): dim_(d) {lineages_ = {};};
    void record(LineageLabel);
    void save(string);
};

inline void LineageRec::record(LineageLabel labels) {
    assert(labels.dim() == dim_);
    lineages_.push_back(labels.labels_);
}

inline void LineageRec::save(string id) {
    // Saver
    std::ofstream writing_file;
    std::string filename = "./results/labels/" + id + ".csv";
    writing_file.open(filename, std::ios::out);

    // Main
    for (vector<int> a : lineages_) { for (int l=0; l<dim_; l++) {
        if (l == dim_-1) {
            writing_file << a[l] << std::endl;
        } else {
            writing_file << a[l] << ",";
        }
    } }
    writing_file.close();
}


namespace PDS {

vector<Lives> zero_life_matrix_generator(int C, int L) {
    vector<Lives> matrix(C);
    for (int c=0; c<C; c++) { matrix[c] = Lives(L, 0.0); }

    return matrix;
}

inline double sum(Lives lives, int L) {
    double sum = 0.0;
    for (int l=0; l<L; l++) { sum += lives[l]; }

    return sum;
}

vector<double> sum_2CLives_d(Lives v1, Lives v2, int dim) {
    vector<double> sum(dim, 0.0);
    for (int i=0; i<dim; i++) { sum[i] = v1[i] + v2[i]; }

    return sum;
}

vector<int> sum_2vecs_int(vector<double> v1, vector<double> v2, int dim) {
    vector<int> sum(dim, 0.0);
    for (int i=0; i<dim; i++) { sum[i] = (int)(v1[i] + v2[i]); }

    return sum;
}

vector<double> sum_lineage(vector<Lives> lives, int C, int L) {
    vector<double> sum(L, 0.0);
    for (int c=0; c<C; c++) { for (int l=0; l<L; l++) {
        sum[l] += lives[c][l];
    } }

    return sum;
}

vector<double> sum_compartment(vector<Lives> lives, int C, int L) {
    vector<double> sum(C, 0.0);
    for (int c=0; c<C; c++) { for (int l=0; l<L; l++) {
        sum[c] += lives[c][l];
    } }

    return sum;
}

Act template_act_generator(int size) {
    Act act(size);
    for (int i=0; i<size; i++) { vector<double> a(size); act[i] = a; }
    return act;
}

Act zero_act_generator(int size) {
    Act act(size);
    for (int i=0; i<size; i++) { vector<double> a(size, 0); act[i] = a; }
    return act;
}

void zero_act(Act& act, int index, int L) {
    for (int l; l<L; l++) { act[l][index] = 0.0; act[index][l] = 0.0; }
}

}


class Env {
public:
    int R, C, FD, L; double S, N; 
    Lives init_lives;
    vector<double> khs, kps;
    LineageLabel init_labels;
    Act init_act;
    double mhh, mpp, mhp;
public:
    Env() {};
    Env(
        int r, int c, double s, int fd, double n, Lives ilives, int l,
        vector<double> kh, vector<double> kp, LineageLabel ilabels, Act iact,
        double mhh_, double mpp_, double mhp_
    ):
        R(r), C(c), S(s), FD(fd), N(n), init_lives(ilives), L(l),
        khs(kh), kps(kp), init_labels(ilabels), init_act(iact),
        mhh(mhh_), mpp(mpp_), mhp(mhp_)
    { assert(L == init_lives.size()); }; // Env's constructor
}; // class Env


namespace PDS {

void display_lives(vector<Lives> lives, int C, int L) {
    vector<double> sum_line = PDS::sum_lineage(lives, C, L);
    std::cout << "Lineage Conc: { ";
    for (int l=0; l<L; l++) {
        std::cout << sum_line[l] << " ";
    }
    std::cout << "}" << std::endl;
}

void display_act(Act act, int L) {
    std::cout << "Act:" << std::endl;
    for (int l1=0; l1<L; l1++) {
    for (int l2=0; l2<L; l2++) { 
        std::cout << act[l1][l2] << " ";
    }
        std::cout << std::endl;
    } 
    std::cout << "}" << std::endl;
}

void display_label(LineageLabel labels, int L) {
    std::cout << "Lineage: { ";
    for (int l=0; l<L; l++) {
        std::cout << labels[l] << " ";
    }
    std::cout << "}" << std::endl;
}

}


#endif  // TYPES_HPP_
