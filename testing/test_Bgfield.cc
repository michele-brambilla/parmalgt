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

class AbeianBgfTest : public CommonBgfTest, public ::testing::Test {
public:
  AbeianBgfTest() {
    for (int i = 0; i < 3; ++i){
      su3B.whr[i*4] = Cplx(r.Rand(), r.Rand());
      bgfB[i] = su3B.whr[i*4];
    }
  };
  SU3 su3B;
  bgf::AbelianBgf bgfB;
};

class TrivialBgfTest : public CommonBgfTest, public ::testing::Test { 
public:
  bgf::TrivialBgf One;
};

TEST(AbelianBgfConstructorTest, ThrowsOnWrongSize){
  std::vector<Cplx> v(4);
  ASSERT_THROW({bgf::AbelianBgf b(v);}, bgf::InitializedWithWrongSize);
}

TEST_F(AbeianBgfTest, ApplyFromRight){
  ASSERT_TRUE( SU3Cmp(bgfB.ApplyFromLeft(A), su3B*A)() );
}

TEST_F(AbeianBgfTest, ApplyFromLeft){
  ASSERT_TRUE( SU3Cmp(bgfB.ApplyFromRight(A), A*su3B)() );
}

TEST_F(TrivialBgfTest, ApplyFromRight){
  ASSERT_TRUE( SU3Cmp(One.ApplyFromRight(A), A)() );
}

TEST_F(TrivialBgfTest, ApplyFromLeft){
  ASSERT_TRUE( SU3Cmp(One.ApplyFromLeft(A), A)() );
}


