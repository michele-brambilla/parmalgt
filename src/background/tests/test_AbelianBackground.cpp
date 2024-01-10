#include <background/AbelianBackground.hpp>
#include <background/TrivialBackground.hpp>
#include "gtest/gtest.h"
#include <common/Helper.hpp>
#include <cstdlib> // for rand

// We need a test for the AbelianBgf class...

// With a nice base class ...

ranlxd::Rand r(23797);

template <class R> static SU3 makeRandomSU3(R& Rand){
    SU3 result;
    using cplx = typename SU3::data_t;
    static double twopi = std::atan(1.0) * 8.0;
    static double soneo3 = std::sqrt( 1./ 3);
    std::array<double, 8> g;
    double rho;
    double theta;
    // get flat distribution
    r.ranlxd(g.begin(), g.end());
    // make gaussian
    for (int i = 0; i < 8; i+=2){
      theta = twopi*g[i];
      rho = std::sqrt( -std::log((1. - g[i+1])) );
      g[i]   = rho*std::cos(theta);
      g[i+1] = rho*std::sin(theta);
    }
    g[7] *= soneo3;

    result[0] = cplx(0,g[7]+g[6]);
    result[1] = cplx(g[1],g[0]);
    result[2] = cplx(g[3],g[2]);
    result[3] = cplx(-g[1],g[0]);
    result[4] = cplx(0,g[7]-g[6]);
    result[5] = cplx(g[5],g[4]);
    result[6] = cplx(-g[3],g[2]);
    result[7] = cplx(-g[5],g[4]);
    result[8] = cplx(0,-g[7]*2);
    return result;
  }

class CommonBgfTest {
public:
  SU3 A{makeRandomSU3(r)};
};

class AbelianBgfTest : public CommonBgfTest, public ::testing::Test {
public:
  AbelianBgfTest() {
    std::array<double, 3*4> random;
    r.ranlxd(random.begin(), random.end());
    for (int i = 0; i < 3; ++i){
      su3B[i*4] = Cplx(random[4*i], random[4*i+1]);
      bgfB[i] = su3B[i*4];
      su3C[i*4] = Cplx(random[4*i+2], random[4*i+3]);
      bgfC[i] = su3C[i*4];
      su3One[i*4] = Cplx(1., 0.);
    }
  };
  SU3 su3B, su3C, su3One;
  bgf::AbelianBgf bgfB, bgfC;
private:

};

class TrivialBgfTest : public CommonBgfTest, public ::testing::Test { 
public:
  bgf::TrivialBgf<SU3, CVector> One;
};


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

TEST_F(AbelianBgfTest, VectorProduct){
  CVector b, c, Bb, cC;
  std::array<double, 3*4> random;
  r.ranlxd(random.begin(), random.end());
  for (int i = 0; i < 3; ++i){
    c[i] =  Cplx(random[4*i], random[4*i+1]);
    b[i] =  Cplx(random[4*i+2], random[4*i+3]);
  }
  cC = bgfC.ApplyFromRight(c);
  Bb = bgfB.ApplyFromLeft(b);
  for (int i = 0; i < 3; ++i){
    ASSERT_DOUBLE_EQ(cC[i].real(), (bgfC[i]*c[i]).real());
    ASSERT_DOUBLE_EQ(cC[i].imag(), (bgfC[i]*c[i]).imag());
    ASSERT_DOUBLE_EQ(Bb[i].real(), (bgfB[i]*b[i]).real());
    ASSERT_DOUBLE_EQ(Bb[i].imag(), (bgfB[i]*b[i]).imag());

  }
}

TEST_F(AbelianBgfTest, CplxProduct){
  std::array<double, 2> random;
  r.ranlxd(random.begin(), random.end());
  Cplx alpha;(random[0],random[1]);
  SU3 D = (bgfB*alpha).ApplyFromLeft(su3One);
  ASSERT_TRUE( SU3Cmp(D, su3B*alpha)() );
}

TEST_F(TrivialBgfTest, ApplyFromRightOnSU3){
  ASSERT_TRUE( SU3Cmp(One.ApplyFromRight(A), A)() );
}

TEST_F(TrivialBgfTest, ApplyFromLeftOnSU3){
  ASSERT_TRUE( SU3Cmp(One.ApplyFromLeft(A), A)() );
}

TEST(AbelianBgf, Arithmetic){
  bgf::AbelianBgf A, B;
  A[0] = 1.2;
  A[1] = 2.3;
  A[2] = 3.4;

  B[0] = 0.2;
  B[1] = 1.3;
  B[2] = 2.4;
  bgf::AbelianBgf C = A+B, D = A-B;
  ASSERT_TRUE(Cmp(C[0], 1.4)());
  ASSERT_TRUE(Cmp(C[1], 3.6)());
  ASSERT_TRUE(Cmp(C[2], 5.8)());
  ASSERT_TRUE(Cmp(D[0], 1.)());
  ASSERT_TRUE(Cmp(D[1], 1.)());
  ASSERT_TRUE(Cmp(D[2], 1.)());
}

TEST(AbelianBgf, KnownValues){
  //  This tests the class
  //     bgf::AbelianBgfFactory
  //  using the formula
  //    V(x0) = V(x0 = 0)* exp{ i a E x0 }
  std::srand(123); 
  // Unit 3x3 matrix
  SU3 su3One;
  su3One[0] = 1;
  su3One[4] = 1;
  su3One[8] = 1;
  // pi / 3
  double pio3 = std::atan(1.)*4./3;
  // do 1000 checks
  for (int _n = 0; _n < 1000; ++_n){
    // random T, L (< 100), eta, nu
    int L = abs((std::rand() * 2)) % 60 + 4;
    int T = L;
    //int T = std::rand() % 100 + 10;
    // initialize the background field
    bgf::AbelianBgfFactory factory(T, L);
    // random x0
    int x0 = std::rand()*T;
    // now calculate 
    //       exp (i E t) = dV
    //           = exp ( - i gamma * x0 * diag(2,-1,-1) )
    double gamma = 1./L/T * pio3;
    SU3 su3dV;
    su3dV[0] = exp(Cplx(0, -2.*gamma*x0));
    su3dV[4] = exp(Cplx(0, gamma*x0));
    su3dV[8] = exp(Cplx(0, gamma*x0));
    // get V(x0) and V(0)
    bgf::AbelianBgf V(factory.get(x0));
    bgf::AbelianBgf V0(factory.get(0));
    EXPECT_TRUE( SU3Cmp(V0.ApplyFromLeft(su3dV),
                        V.ApplyFromLeft(su3One))());
  }
}

// Test tree level formula for E,

TEST(AbelianBgf, E){
  int L = 8;
  int T = 8;
  // initialize
  bgf::get_abelian_bgf(0,0,T,L,0);
  SU3 C;
  C(0,0) = Cplx( 0, 1./ L);
  C(1,1) = Cplx( 0, -.5/ L);
  C(2,2) = Cplx( 0, -.5/ L);
  Cplx E = 0;
  for (int k = 1; k < 4; ++k){
    E -= (bgf::get_abelian_bgf(0, k)*
          bgf::get_abelian_bgf(0, 0)*
          bgf::get_abelian_bgf(1, k).dag()*
          bgf::get_abelian_bgf(0, 0).dag()).ApplyFromRight(C).tr();
    E -= (bgf::get_abelian_bgf(T, k).dag()*
          bgf::get_abelian_bgf(T-1, 0).dag()*
          bgf::get_abelian_bgf(T-1, k)*
          bgf::get_abelian_bgf(T-1, 0).dag()).ApplyFromRight(C).tr();
  }
  E *= L*L*L*2;
  double pio3 = std::atan(1.)*4./3;
  double k = 12*L*L*(std::sin(pio3/L/L) + std::sin(pio3*2/L/L));
  EXPECT_DOUBLE_EQ( E.real(), k);


  L = 8;
  int s = -1;
  T = L - s;
  C(0,0) = Cplx( 0, 1./ L);
  C(1,1) = Cplx( 0, -.5/ L);
  C(2,2) = Cplx( 0, -.5/ L);
  bgf::AbelianBgfFactory factory(T, L, s);
  E = 0;
  E -= (bgf::AbelianBgf(factory.get(0))*
        bgf::AbelianBgf(factory.get(1)).dag()).ApplyFromRight(C).tr();
  E -= (bgf::AbelianBgf(factory.get(T)).dag()*
        bgf::AbelianBgf(factory.get(T-1))).ApplyFromRight(C).tr();
  E *= L*L*L*2*3;
  EXPECT_DOUBLE_EQ( E.real(), 37.6945384607827);
  
  s = 1;
  T = L - s;
  C(0,0) = Cplx( 0, 1./ L);
  C(1,1) = Cplx( 0, -.5/ L);
  C(2,2) = Cplx( 0, -.5/ L);
  factory = bgf::AbelianBgfFactory(T, L, s);
  E = 0;
  E -= (bgf::AbelianBgf(factory.get(0))*
        bgf::AbelianBgf(factory.get(1)).dag()).ApplyFromRight(C).tr();
  E -= (bgf::AbelianBgf(factory.get(T)).dag()*
        bgf::AbelianBgf(factory.get(T-1))).ApplyFromRight(C).tr();
  E *= L*L*L*2*3;
  ASSERT_DOUBLE_EQ( E.real(), 37.69169953984173);
}
  

