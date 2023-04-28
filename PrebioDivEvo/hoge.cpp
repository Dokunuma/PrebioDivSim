#if __INTELLISENSE__
#undef __ARM_NEON
#undef __ARM_NEON__
#endif

#include <iostream>
#include <typeinfo>
#include <fstream>
#include <Eigen/Dense>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::VectorXi;

int main()
{
    VectorXd v1 = VectorXd::Random(3);
    std::cout << v1 << std::endl;
    VectorXd v2 = VectorXd::Random(3);
    std::cout << v2 << std::endl;
    VectorXi v3 = (v1 + v2).cast<int>();
    std::cout << v3 << std::endl;
}