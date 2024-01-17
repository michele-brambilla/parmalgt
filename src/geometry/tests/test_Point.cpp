#include <geometry/Point.hpp>
#include "gtest/gtest.h"

const int DIM = 4;

using namespace pt;

template <int D>
bool test_is_good_true_if_0_le_direction() {
    Direction<D> dir;
    while (dir.is_good()) {
        --dir;
    }
    return int(dir) == -1;
}

template <int D>
bool test_is_good_true_if_direction_lt_dim() {
    Direction<D> dir;
    while (dir.is_good()) {
        ++dir;
    }
    return int(dir) == D;
}

TEST(Direction, test_is_good_true_if_0_le_direction_lt_DIM) {
    EXPECT_TRUE(test_is_good_true_if_0_le_direction<1>());
    EXPECT_TRUE(test_is_good_true_if_direction_lt_dim<1>());

    EXPECT_TRUE(test_is_good_true_if_0_le_direction<2>());
    EXPECT_TRUE(test_is_good_true_if_direction_lt_dim<2>());

    EXPECT_TRUE(test_is_good_true_if_0_le_direction<3>());
    EXPECT_TRUE(test_is_good_true_if_direction_lt_dim<3>());

    EXPECT_TRUE(test_is_good_true_if_0_le_direction<4>());
    EXPECT_TRUE(test_is_good_true_if_direction_lt_dim<4>());

    EXPECT_TRUE(test_is_good_true_if_0_le_direction<5>());
    EXPECT_TRUE(test_is_good_true_if_direction_lt_dim<5>());

    EXPECT_TRUE(test_is_good_true_if_0_le_direction<10>());
    EXPECT_TRUE(test_is_good_true_if_direction_lt_dim<10>());

    EXPECT_TRUE(test_is_good_true_if_0_le_direction<100>());
    EXPECT_TRUE(test_is_good_true_if_direction_lt_dim<100>());
}