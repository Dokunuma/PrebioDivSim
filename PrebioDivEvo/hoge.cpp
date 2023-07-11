#include "src/Types.hpp"

int main()
{
    vector<double> init_lives{100.0, 0.0, 0.0};
    Live lives(3000, 10);

    lives.filled(init_lives);
    std::cout << lives.values(0) << std::endl;
}