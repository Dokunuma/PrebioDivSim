#ifndef SAMPLE_HPP_
#define SAMPLE_HPP_

#include <random>
#include <vector>
#include <iostream>
#include <cstdint>
#include <algorithm>
#include <stdexcept>
#include <unordered_set>
#include <functional>

namespace PDS {

std::mt19937 create_rand_engine(){
    std::random_device rnd;
    std::vector<std::uint_least32_t> v(10);// 初期化用ベクタ
    std::generate(v.begin(), v.end(), std::ref(rnd));// ベクタの初期化
    std::seed_seq seed(v.begin(), v.end());
    return std::mt19937(seed);// 乱数エンジン
}

std::vector<int> make_rand_array(const size_t size, int rand_min, int rand_max) {
    if (rand_min > rand_max) std::swap(rand_min, rand_max);
    const size_t max_min_diff = static_cast<size_t>(rand_max - rand_min + 1);
    if (max_min_diff < size) throw std::runtime_error("引数が異常です");

    std::unordered_set<int> tmp;
    auto engine = create_rand_engine();
    std::uniform_int_distribution<int> distribution(rand_min, rand_max);
    while (tmp.size() < size) tmp.insert(distribution(engine));
    return std::vector<int>(tmp.begin(), tmp.end());
}

std::vector<int> make_rand_array(
    const size_t size, int rand_min, int rand_max, std::mt19937& engine
) {
    if (rand_min > rand_max) std::swap(rand_min, rand_max);
    const size_t max_min_diff = static_cast<size_t>(rand_max - rand_min + 1);
    if (max_min_diff < size) throw std::runtime_error("引数が異常です");

    std::unordered_set<int> tmp;
    std::uniform_int_distribution<int> distribution(rand_min, rand_max);
    while (tmp.size() < size) tmp.insert(distribution(engine));
    
    return std::vector<int>(tmp.begin(), tmp.end());
}

}

#endif  // UTILS_HPP_
