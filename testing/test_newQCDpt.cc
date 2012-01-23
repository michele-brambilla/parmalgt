#include "newQCDpt.h"
#include <Background.h>
#include "gtest/gtest.h"
#include "Helper.h"
#include <cstdlib> // for rand

const int PT_ORD = 8;
const int AL_ORD = 10;

typedef BGptSU3<bgf::AbelianBgf, AL_ORD, PT_ORD> ptSU3;
typedef BGptCVector<AL_ORD, PT_ORD> ptCVector;
typedef BGptGluon<bgf::AbelianBgf, AL_ORD, PT_ORD, 4> ptGluon;
typedef BGptSpinColor<AL_ORD, PT_ORD, 4> ptSpinColor;

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
    MyPtSU3A(), 
    MyPtSU3B() 
  { 
    MyPtSU3A.randomize();
    MyPtSU3B.randomize();
  }
  ptSU3 MyPtSU3A, MyPtSU3B;
};

// Multiplication of two ptSU3s of the form
// A = V(t1) + g_0^2 A1 + ... ,
// B = V(t2) + g_0^2 B1 + ... , and
// One = 1

TEST_F(AbelianBgfTest, SimpleMultiply){
  ptSU3 One(bgf::unit()); // unit matrix
  // to be extra safe, check multiplication from left and right
  ptSU3 ACopy = One*MyPtSU3A;
  ptSU3 BCopy = MyPtSU3B*One;
  for (int i = 0; i < PT_ORD; ++i){
    ASSERT_TRUE( SU3Cmp(ACopy[i], MyPtSU3A[i])() );
    ASSERT_TRUE( SU3Cmp(BCopy[i], MyPtSU3B[i])() );
  }
}

TEST(AbelianBgf, Multiply){
  ptSU3 A(bgf::unit()), B(bgf::unit());
  SU3 One, Zero;
  for (int i = 0; i < 3; ++i)
    One(i,i) = Cplx(1,0);
  A[2] = One*2;
  B[4] = One*3;
  ptSU3 AB = A * B;
  ASSERT_TRUE( SU3Cmp(AB[0], Zero)() );
  ASSERT_TRUE( SU3Cmp(AB[1], Zero)() );
  ASSERT_TRUE( SU3Cmp(AB[2], One*2)() );
  ASSERT_TRUE( SU3Cmp(AB[3], Zero)() );
  ASSERT_TRUE( SU3Cmp(AB[4], One*3)() );
  ASSERT_TRUE( SU3Cmp(AB[5], Zero)() );
  ASSERT_TRUE( SU3Cmp(AB[6], Zero)() );
  ASSERT_TRUE( SU3Cmp(AB[7], One*6)() );
}

// Test the addition,
// A+A == 2*A

TEST_F(AbelianBgfTest, Add){
  ptSU3 A_times_2 = MyPtSU3A*2.;
  ptSU3 A_plus_A = MyPtSU3A + MyPtSU3A;
  for (int i = 0; i < PT_ORD; ++i)
    ASSERT_TRUE( SU3Cmp(A_times_2[i], A_plus_A[i])() );
  ASSERT_TRUE( A_times_2.bgf() ==  A_plus_A.bgf() );
}


// Testing the scalar multiplicatiohn

TEST_F(AbelianBgfTest, ScalarMultiply){
  Cplx alpha(r.Rand(), r.Rand());
  Cplx beta(r.Rand(), r.Rand());
  ptSU3 beta_B = MyPtSU3B;
  beta_B *= beta;
  for (int i = 0; i < PT_ORD; ++i){
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

TEST(BGptCVector, ArithmeticPlusSelf){
  ptCVector u, v, w;
  for( ptCVector::iterator i = v.begin(), j = w.begin();
         i != v.pt_end(); ++i, ++j)
    for (int k = 0; k < 3; ++k){
      i->whr[k] = Cplx(r.Rand(), r.Rand());
      j->whr[k] = i->whr[k]*3;
    }
  u = v+v;
  u += v;
  for( ptCVector::iterator i = u.begin(), j = w.begin();
         i != u.pt_end(); ++i, ++j)
    ASSERT_TRUE(*i == *j);
}

TEST(BGptCVector, ProductWithSU3){
  ptCVector v, w;
  SU3 One;
  CVector OneV, ZeroV;
  ptSU3 A(bgf::unit());
  for (int i = 0; i < 3; ++i){
    One(i,i) = Cplx(1,0);
    OneV.whr[i] = Cplx(1,0);
  }
  A[3] = One*2;
  for (int i = 0; i < 3; ++i)
    v[1].whr[i] = Cplx(3,0);
  w = v*A;
  ASSERT_TRUE(w[0] == ZeroV);
  ASSERT_TRUE(w[1] == OneV*Cplx(3,0));
  ASSERT_TRUE(w[2] == ZeroV);
  ASSERT_TRUE(w[3] == ZeroV);
  ASSERT_TRUE(w[4] == ZeroV);
  ASSERT_TRUE(w[5] == OneV*Cplx(6,0));
  ASSERT_TRUE(w[6] == ZeroV);
}

TEST(BGptSpinColor, ProductWithGluon){
  ptGluon U;
  SU3 One;
  CVector OneV, ZeroV;
  ptSU3 A(bgf::unit());
  ptSpinColor psi, chi;
  for (int i = 0; i < 3; ++i){
    One(i,i) = Cplx(1,0);
    OneV.whr[i] = Cplx(1,0);
  }
  for (int i = 0; i < 4; ++i){
    U[i].bgf() = bgf::unit();
    U[i][i] = One*i;
    psi[i][5] = OneV*i;
  }
  chi = psi*U; // <-- check this
  for (int i = 0; i < 4; ++i)
    for (int j = 0; j <= PT_ORD; ++j)
      if (i && j == 5)
        ASSERT_TRUE(chi[i][j] == OneV*Cplx(i));
      else if (i && j == 5 + i + 1)
        ASSERT_TRUE(chi[i][j] == OneV*Cplx(i)*Cplx(i));
      else
        ASSERT_TRUE(chi[i][j] == ZeroV);
}

// we need a main function here because we have to initialize the
// background field class...

int main(int argc, char **argv) {
  bgf::get_abelian_bgf(0,0,T, L, eta, nu);
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
