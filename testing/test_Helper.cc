#include "Helper.h"
#include "gtest/gtest.h"

MyRand r(32490714);

TEST(SU3CmpTest, EqualMatrices){
  SU3 A,B;
  for (int i = 0; i < 9; ++i){
    A.whr[i] = Cplx(r.Rand(), r.Rand());
    B.whr[i] = A.whr[i];
  }
  ASSERT_TRUE(SU3Cmp(A,B)(0.0));
}

TEST(SU3CmpTest, NotEqualMatrices){
  SU3 A,B;
  for (int i = 0; i < 9; ++i){
    A.whr[i] = Cplx(r.Rand(), r.Rand());
    B.whr[i] = A.whr[i];
  }
  for (int i = 8; i >= 0; --i){
    A = B;
    A.whr[i] = Cplx(r.Rand(), r.Rand());
    ASSERT_FALSE(SU3Cmp(A,B)());
  }
}
