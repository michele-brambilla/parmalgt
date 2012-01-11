int PTORD = 6;

#include "newQCDpt.h"
#include <Background.h>
#include "gtest/gtest.h"
#include "Helper.h"
#include <cstdlib> // for rand

MyRand r(23797);

int L = std::rand() % 100;
int T = std::rand() % 100;
double eta = r.Rand();
double nu = r.Rand();
// initialize the background field


// Simple test base class defining two random 
// perturbative SU3s ,,,

class AbelianBgfTest: public ::testing::Test {
public:
  AbelianBgfTest() : 
    MyPtSU3A(bgf::AbelianBgf(T*r.Rand())), 
    MyPtSU3B(bgf::AbelianBgf(T*r.Rand())) 
  { 
    MyPtSU3A.randomize();
    MyPtSU3B.randomize();
  }
  BGptSU3<bgf::AbelianBgf> MyPtSU3A, MyPtSU3B;
};

// Multiplication of two ptSU3s of the form
// A = V(t1) + g_0^2 A1 + ... ,
// B = V(t2) + g_0^2 B1 + ... , and
// One = 1

TEST_F(AbelianBgfTest, SimpleMultiply){
  std::vector<Cplx> v(3, 1);
  bgf::AbelianBgf UnitBgf(v); // unit bgf diag(1,1,1)
  BGptSU3<bgf::AbelianBgf> One(UnitBgf); // unit matrix
  // to be extra safe, check multiplication from left and right
  BGptSU3<bgf::AbelianBgf> ACopy = One*MyPtSU3A;
  BGptSU3<bgf::AbelianBgf> BCopy = MyPtSU3B*One;
  for (int i = 0; i < PTORD; ++i){
    ASSERT_TRUE( SU3Cmp(ACopy[i], MyPtSU3A[i])() );
    ASSERT_TRUE( SU3Cmp(BCopy[i], MyPtSU3B[i])() );
  }
}

// Testing the scalar multiplicatiohn

TEST_F(AbelianBgfTest, ScalarMultiply){
  Cplx alpha(r.Rand(), r.Rand());
  Cplx beta(r.Rand(), r.Rand());
  BGptSU3<bgf::AbelianBgf> beta_B = MyPtSU3B;
  beta_B *= beta;
  for (int i = 0; i < PTORD; ++i){
    SU3 a = MyPtSU3A[i];
    SU3 b = MyPtSU3B[i];
    // multipy 'by hand'
    for (SU3::iterator ait = a.begin(), bit = b.begin();
         ait != a.end(); ++ait, ++bit){
      *ait *= alpha;
      *bit *= beta;
    }
    EXPECT_TRUE ( SU3Cmp( (MyPtSU3A*alpha)[i], a)());
    EXPECT_TRUE ( SU3Cmp( beta_B[i], b)());
  }
}

// we need a main function here because we have to initialize the
// background field class...

int main(int argc, char **argv) {
  bgf::AbelianBgf::init(T, L, eta, nu);
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
