#pragma once
#include <random>

namespace rng_utils
{
    extern const int SEED = 1998; // Shoutout!
    static std::random_device rd;
    static std::mt19937 rng(SEED);
    static std::uniform_real_distribution<double> dis(0.0, 1.0);
}