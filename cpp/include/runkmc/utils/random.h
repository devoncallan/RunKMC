#pragma once
#include <random>

namespace rng
{

    namespace core
    {
        constexpr int SEED = 1998; // Shoutout!
        static std::random_device rd;
        static std::mt19937_64 generator(SEED);
        static std::uniform_real_distribution<double> dis(0.0, 1.0);
    }

    // Returns a random double in (0, 1]
    static inline double rand()
    {
        return core::dis(core::generator) + 1e-40;
    }

    // Returns a random index in [0, max)
    static inline size_t randIndex(size_t max)
    {
        return static_cast<size_t>(core::dis(core::generator) * max);
    }

    // Returns a random index sampled from a discrete weighted distribution
    template <typename T>
    static inline size_t randIndexWeighted(const std::vector<T> &weights)
    {
        std::discrete_distribution<size_t> discrete_dis(weights.begin(), weights.end());
        return discrete_dis(core::generator);
    }
}