#ifndef UTILS_HPP_
#define UTILS_HPP_

#include "Types.hpp"
#include "Sample.hpp"
#include <sys/time.h>
#include <ctime>
#include <chrono>
#include <random>
#include <algorithm>

using std::chrono::duration_cast;
using std::chrono::milliseconds;
using std::chrono::system_clock;

namespace PDS {

vector<int> range(int size) {
    vector<int> vec(size, 0);

    for (int i=0; i<size; i++) {
        vec[i] = i;
    }

    return vec;
}

vector<int> selection_mapper(int C, int size, std::mt19937& engine) {
    vector<int> selected = PDS::make_rand_array(size, 0, C-1, engine);
    return selected;
}

vector<vector<int>> fusdiv_mapper(int C, int size, std::mt19937& engine) {
    vector<vector<int>> fusdived = {};
    for (int s=0; s<size; s++) {
        vector<int> fd = PDS::make_rand_array(2, 0, C-1, engine);
        fusdived.push_back(fd);
    }
    
    return fusdived;
}

string get_time_str() {
    string mtime;
    struct timespec ts;
    struct tm t;

    // Get epoch time
    clock_gettime(CLOCK_REALTIME, &ts);
    // Convert into local and parsed time
    localtime_r(&ts.tv_sec, &t);
    // Create string with strftime
    char buf[32];
    strftime(buf, 32, "%Y-%m-%d-%H-%M-%S", &t);
    // Add milli-seconds with snprintf
    char output[32];
    const int msec = ts.tv_nsec / 1000000;
    snprintf(output, 32, "%s-%03d", buf, msec);
    // Result
    mtime = string(output);

    return mtime;
}

string random_code_generator(int length) {
    string code;
    std::mt19937 mt{std::random_device{}()};
    std::uniform_int_distribution<int> dist(0, 25); // 26 letters of a~z

    for (int i = 0; i < length; i++) { code += char(dist(mt) + 'a'); }

    return code;
}


// DEBUG functions
void debug_selection(int R, vector<vector<int>> killed) {
    for (int rd=0; rd<R; rd++) {
        if (rd==0) {
            std::cout << "### Index for Selection ###"<< std::endl;
        }
        std::cout << "Round " << rd+1 << "(" << killed[rd].size() << ") {";
        for (int kill : killed[rd]) { std::cout << kill << " "; }
        std::cout << "}" << std::endl;
}  } // funcdebug_selection

void debug_fusion_division(int R, vector<vector<vector<int>>> fusdived) {
    for (int rd=0; rd<R; rd++) {
        if (rd==0) {
            std::cout << "### Index for Fusion-Division ###"<< std::endl;
        }
        std::cout << "Round " << rd+1 << "(" << fusdived[rd].size() << ") ";
        for (vector<int> fd : fusdived[rd]) {
            std::cout << "{" << fd[0] << " " << fd[1] << "}";
    } std::cout << std::endl; 
}  } // func debug_fusion_division


} // namespce PDS

#endif  // UTILS_HPP_
