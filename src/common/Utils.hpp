# ifndef _COMMON_UTILS_HPP
# define _COMMON_UTILS_HPP

#include <types/Types.hpp>

#include <random>

static std::random_device rd;
static std::mt19937 gen{rd()};
static const double soneo3 = std::sqrt(1. / 3);

SU3 makeRandomSU3() {
    SU3 result;
    using cplx = typename SU3::data_t;
    std::array<double, 8> g;
    std::normal_distribution<double> distrib(0, 1);
    for (auto &item : g) {
        item = distrib(gen);
    }
    g[7] *= soneo3;

    result[0] = cplx(0, g[7] + g[6]);
    result[1] = cplx(g[1], g[0]);
    result[2] = cplx(g[3], g[2]);
    result[3] = cplx(-g[1], g[0]);
    result[4] = cplx(0, g[7] - g[6]);
    result[5] = cplx(g[5], g[4]);
    result[6] = cplx(-g[3], g[2]);
    result[7] = cplx(-g[5], g[4]);
    result[8] = cplx(0, -g[7] * 2);
    return result;
}



#endif