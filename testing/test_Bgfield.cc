#include "Background.h"
#include "gtest/gtest.h"
#include "Helper.h"

// We need a test for the AbelianBgf class...

// With a nice base class ...

MyRand r(23797);




class CommonBgfTest {
public:
  CommonBgfTest() {
    for (int i = 0; i < 9; ++i)
      A.whr[i] = Cplx(r.Rand(), r.Rand());
  }
  SU3 A;
};

class AbeianBgfUnitBTest : public CommonBgfTest, public ::testing::Test {
public:
  AbeianBgfUnitBTest() {
    for (int i = 0; i < 3; ++i){
      su3B.whr[i*4] = Cplx(r.Rand(), r.Rand());
      bgfB[i] = su3B.whr[i*4];
    }
  };
  SU3 su3B;
  bgf::AbelianBgf bgfB;
};

TEST_F(AbeianBgfUnitBTest, ApplyFromRight){
  ASSERT_TRUE( SU3Cmp(bgfB.ApplyFromLeft(A), su3B*A)() );
};

TEST_F(AbeianBgfUnitBTest, ApplyFromLeft){
  ASSERT_TRUE( SU3Cmp(bgfB.ApplyFromRight(A), A*su3B)() );
};

