#include "Background.h"
#include "gtest/gtest.h"
#include "Helper.h"
#include <cstdlib> // for rand

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

class AbelianBgfTest : public CommonBgfTest, public ::testing::Test {
public:
  AbelianBgfTest() {
    for (int i = 0; i < 3; ++i){
      su3B.whr[i*4] = Cplx(r.Rand(), r.Rand());
      bgfB[i] = su3B.whr[i*4];
      su3C.whr[i*4] = Cplx(r.Rand(), r.Rand());
      bgfC[i] = su3C.whr[i*4];
      su3One.whr[i*4] = Cplx(1., 0.);
    }
  };
  SU3 su3B, su3C, su3One;
  bgf::AbelianBgf bgfB, bgfC;
};

class TrivialBgfTest : public CommonBgfTest, public ::testing::Test { 
public:
  bgf::TrivialBgf One;
};

TEST(AbelianBgfConstructorTest, ThrowsOnWrongSize){
  std::vector<Cplx> v(4);
  ASSERT_THROW({bgf::AbelianBgf b(v);}, bgf::InitializedWithWrongSize);
}
TEST_F(AbelianBgfTest, CopyConstructor){
  bgf::AbelianBgf Bcopy(bgfB);
  bgf::AbelianBgf Ccopy = bgfC;
  ASSERT_TRUE(Bcopy == bgfB);
  ASSERT_TRUE(Ccopy == bgfC);
  ASSERT_FALSE(Bcopy == Ccopy);
}
TEST_F(AbelianBgfTest, ApplyFromRight){
  ASSERT_TRUE( SU3Cmp(bgfB.ApplyFromLeft(A), su3B*A)() );
}

TEST_F(AbelianBgfTest, ApplyFromLeft){
  ASSERT_TRUE( SU3Cmp(bgfB.ApplyFromRight(A), A*su3B)() );
}

TEST_F(AbelianBgfTest, BgfProduct){
  SU3 D = (bgfB*bgfC).ApplyFromLeft(su3One);
  ASSERT_TRUE( SU3Cmp(D, su3B*su3C)() );
}

TEST_F(AbelianBgfTest, CplxProduct){
  Cplx alpha(r.Rand(), r.Rand());
  SU3 D = (bgfB*alpha).ApplyFromLeft(su3One);
  ASSERT_TRUE( SU3Cmp(D, su3B*alpha)() );
}

TEST_F(TrivialBgfTest, ApplyFromRight){
  ASSERT_TRUE( SU3Cmp(One.ApplyFromRight(A), A)() );
}

TEST_F(TrivialBgfTest, ApplyFromLeft){
  ASSERT_TRUE( SU3Cmp(One.ApplyFromLeft(A), A)() );
}

TEST(AbelianBgf, KnownValues){
  //  This tests the constructor
  //     bgf::AbelianBgf::AbelianBgf(int t),
  //  using the formula
  //    V(x0) = V(x0 = 0)* exp{ i a E x0 }
  std::srand(123); 
  // Unit 3x3 matrix
  SU3 su3One;
  su3One.whr[0] = 1;
  su3One.whr[4] = 1;
  su3One.whr[8] = 1;
  // pi / 3
  double pio3 = std::atan(1.)*4./3;
  // do 1000 checks
  for (int _n = 0; _n < 1000; ++_n){
    // random T, L (< 100), eta, nu
    int L = std::rand() % 100;
    int T = std::rand() % 100;
    double eta = r.Rand();
    double nu = r.Rand();
    // initialize the background field
    bgf::AbelianBgf::init(T, L, eta, nu);
    // random x0
    int x0 = r.Rand()*T;
    // now calculate 
    //       exp (i E t) = dV
    //           = exp ( - i gamma * x0 * diag(2,-1,-1) )
    double gamma = 1./L/T * (eta + pio3);
    SU3 su3dV;
    su3dV.whr[0] = exp(Cplx(0, -2.*gamma*x0));
    su3dV.whr[4] = exp(Cplx(0, gamma*x0));
    su3dV.whr[8] = exp(Cplx(0, gamma*x0));
    // get V(x0) and V(0)
    bgf::AbelianBgf V(x0);
    bgf::AbelianBgf V0(0);
    EXPECT_TRUE( SU3Cmp(V0.ApplyFromLeft(su3dV),
                        V.ApplyFromLeft(su3One))());
  }
}
