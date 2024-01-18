#include "Common.hpp"

#include <random>

#include "gtest/gtest.h"

struct TrivialBgfTest : public CommonBgfTest, public ::testing::Test {
    bgf::TrivialBgf One;
};

TEST_F(TrivialBgfTest, ApplyFromRightOnSU3) {
    ASSERT_TRUE(SU3Cmp(One.ApplyFromRight(A), A)());
}

TEST_F(TrivialBgfTest, ApplyFromLeftOnSU3) {
    ASSERT_TRUE(SU3Cmp(One.ApplyFromLeft(A), A)());
}
