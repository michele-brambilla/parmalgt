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

TEST_F(AbelianBgfTest, SimpleMultiply){
  std::vector<Cplx> v(3, 1);
  bgf::AbelianBgf UnitBgf(v);
  BGptSU3<bgf::AbelianBgf> One(UnitBgf);
  BGptSU3<bgf::AbelianBgf> ACopy = One*MyPtSU3A;
  BGptSU3<bgf::AbelianBgf> BCopy = MyPtSU3B*One;
  for (int i = 0; i < PTORD; ++i){
    ASSERT_TRUE( SU3Cmp(ACopy[i], MyPtSU3A[i])() );
    ASSERT_TRUE( SU3Cmp(BCopy[i], MyPtSU3B[i])() );
  }
}

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

TEST(SU3, Multiplication){
  SU3 A, Acpy, alphaA;
  Cplx alpha(r.Rand(), r.Rand());
  for (SU3::iterator ait = A.begin(), 
         aait = alphaA.begin(),
         acpit = Acpy.begin(); 
       ait != A.end(); ++ait, ++aait, ++acpit){
    *ait = Cplx(r.Rand(), r.Rand());
    *acpit = *ait;
    *aait = *ait*alpha;
  }
  Acpy *= alpha;
  EXPECT_TRUE ( SU3Cmp( alphaA, alpha*A)() );
  EXPECT_TRUE ( SU3Cmp( alphaA, A*alpha)() );
  EXPECT_TRUE ( SU3Cmp( alphaA, Acpy)() );
}

int main(int argc, char **argv) {
  bgf::AbelianBgf::init(T, L, eta, nu);
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
