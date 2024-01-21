#include <qcd/newQCDpt.hpp>
#include <background/AbelianBackground.hpp>
#include <common/Helper.hpp>
#include <common/Utils.hpp>


#include "gtest/gtest.h"


const int MYORD = 10;

using ptSU3 = BGptSU3<bgf::AbelianBgf, MYORD>;
using ptsu3 = ptt::PtMatrix<MYORD>;
//typedef BGptCVector<MYORD> ptCVector;
using ptGluon = BGptGluon<bgf::AbelianBgf, MYORD, 4>;
//typedef BGptSpinColor<MYORD, 4> ptSpinColor;

int L = std::rand() % 100;
int T = std::rand() % 100;
// initialize the background field



// //TEST(BGptCVector, ArithmeticPlusSelf){
// //  ptCVector u, v, w;
// //  for( ptCVector::iterator i = v.begin(), j = w.begin();
// //         i != v.end(); ++i, ++j)
// //    for (int k = 0; k < 3; ++k){
// //      i->whr[k] = Cplx(r.Rand(), r.Rand());
// //      j->whr[k] = i->whr[k]*3;
// //    }
// //  u = v+v;
// //  u += v;
// //  for( ptCVector::iterator i = u.begin(), j = w.begin();
// //         i != u.end(); ++i, ++j)
// //    ASSERT_TRUE(*i == *j);
// //}
// //
// //TEST(BGptCVector, ProductWithSU3){
// //  ptCVector v, w;
// //  SU3 One;
// //  CVector OneV, ZeroV;
// //  ptSU3 A(bgf::unit<bgf::AbelianBgf>());
// //  for (int i = 0; i < 3; ++i){
// //    One(i,i) = Cplx(1,0);
// //    OneV.whr[i] = Cplx(1,0);
// //  }
// //  A[3] = One*2;
// //  for (int i = 0; i < 3; ++i)
// //    v[1].whr[i] = Cplx(3,0);
// //  w = v*A;
// //  ASSERT_TRUE(w[0] == ZeroV);
// //  ASSERT_TRUE(w[1] == OneV*Cplx(3,0));
// //  ASSERT_TRUE(w[2] == ZeroV);
// //  ASSERT_TRUE(w[3] == ZeroV);
// //  ASSERT_TRUE(w[4] == ZeroV);
// //  ASSERT_TRUE(w[5] == OneV*Cplx(6,0));
// //  ASSERT_TRUE(w[6] == ZeroV);
// //}

// //TEST(BGptSpinColor, ProductWithGluon){
// //  ptGluon U;
// //  SU3 One;
// //  CVector OneV, ZeroV;
// //  ptSU3 A(bgf::unit<bgf::AbelianBgf>());
// //  ptSpinColor psi, chi;
// //  for (int i = 0; i < 3; ++i){
// //    One(i,i) = Cplx(1,0);
// //    OneV.whr[i] = Cplx(1,0);
// //  }
// //  for (int i = 0; i < 4; ++i){
// //    U[i].bgf() = bgf::unit<bgf::AbelianBgf>();
// //    U[i][i] = One*i;
// //    psi[i][5] = OneV*i;
// //  }
// //  chi = psi*U; // <-- check this
// //  for (int i = 0; i < 4; ++i)
// //    for (int j = 0; j <= MYORD; ++j)
// //      if (i && j == 5)
// //        ASSERT_TRUE(chi[i][j] == OneV*Cplx(i));
// //      else if (i && j == 5 + i + 1)
// //        ASSERT_TRUE(chi[i][j] == OneV*Cplx(i)*Cplx(i));
// //      else
// //        ASSERT_TRUE(chi[i][j] == ZeroV);
// //}

TEST(BGptSU3Test, dag){
 ptSU3 A{bgf::random()};
 ptSU3 B;
 A.randomize();
 B = dag(A);
 for (int r = 0; r < MYORD; ++r)
   for (int i = 0; i < 3; ++i)
     for (int j = 0; j < 3; ++j){
       ASSERT_DOUBLE_EQ( B[r](i,j).real(), A[r](j,i).real() );
     }
}

// TEST(BGptSU3Test, exp){
//   ptt::PtMatrix<MYORD> A = ptt::get_random_pt_matrix<MYORD>();
//   ptSU3 B = exp<bgf::AbelianBgf, MYORD>(A);
//   // calculate the first few terms by hand ...
//   // tree level
//   ASSERT_TRUE(B.bgf() == bgf::unit<bgf::AbelianBgf>());
//   // one loop
//   ASSERT_TRUE(SU3Cmp(A[0],B[0])());
//   // two loop
//   SU3 e = A[0]*A[0]*0.5 + A[1];  
//   ASSERT_TRUE(SU3Cmp(e,B[1])());
//   // three loop
//   e = A[0]*A[0]*A[0]/6. + A[0]*A[1]*.5 + A[1]*A[0]*.5 + A[2];
//   double eps = std::numeric_limits<double>::epsilon() * 10;
//   ASSERT_TRUE(SU3Cmp(e,B[2])(eps));
// }

// Check if exp(A) * exp(-A) = 1 holds

TEST(BGptSU3Test, expUnit){
  for (int _n = 0; _n < 1000; ++_n){
    ptt::PtMatrix<MYORD> A = ptt::get_random_pt_matrix<MYORD>();
    ptSU3 B = exp<bgf::AbelianBgf, MYORD>(A);
    ptSU3 C = exp<bgf::AbelianBgf, MYORD>(A * -1.);
    ptSU3 D = B*C;
    SU3 zero;
    // calculate the first few terms by hand ...
    // tree level
    ASSERT_TRUE(D.bgf() == bgf::unit<bgf::AbelianBgf>());
    for (int i = 0; i < 7; ++i)
      ASSERT_TRUE( SU3Cmp( zero, D[i] )(2e-14*i) );
  }
}

TEST(BGptSU3Test, inverse){
  ptSU3 Z;
  Z.ptU() = ptt::get_random_pt_matrix<MYORD>();
  ptSU3 W = inverse(Z);
  ptSU3 unit = Z*W;
  ASSERT_TRUE( unit.bgf() == bgf::unit<bgf::AbelianBgf>() );
  for (int i = 0; i < MYORD; ++i)
    ASSERT_TRUE( SU3Cmp( unit[i], SU3() )(1e-10) ) << i;
}

TEST(BGptSU3Test, expUnitMore){
  for (int _n = 0; _n < 1000; ++_n){
    ptt::PtMatrix<MYORD> A = ptt::get_random_pt_matrix<MYORD>();
    ptSU3 a, b;
    a.randomize(); a.bgf() = bgf::random();
    b.randomize(); b.bgf() = bgf::random();
    ptSU3 d = b;
    ptSU3 B = exp<bgf::AbelianBgf, MYORD>(A);
    ptSU3 C = exp<bgf::AbelianBgf, MYORD>(A * -1.);
    ptSU3 D = (B)*(C*b);
    // calculate the first few terms by hand ...
    // tree level
    ASSERT_TRUE(D.bgf() == d.bgf());
    for (int i = 0; i < 7; ++i)
      ASSERT_TRUE( SU3Cmp( d[i], D[i] )(1e-13*(i+1)) ) << i;
  }
}

TEST(BGptSU3Test, Multiplication){
  for (int _n = 0; _n < 1000; ++_n){
    ptSU3 a, b, c, d, e;
    a.randomize(); a.bgf() = bgf::random();
    b.randomize(); b.bgf() = bgf::random();
    c.randomize(); c.bgf() = bgf::random();
    d = (a*b)*c; e = a*(b*c);
    for (int i = 0; i < MYORD; ++i)
      ASSERT_TRUE( SU3Cmp( d[i], e[i] )(1e-13*(i+1)) ) << i;
  }
}    

TEST(BGptSU3Test, logThrows){
  ptSU3 A(bgf::random());
  if (IsZero<do_debug, bgf::AbelianBgf>::debug_on)
    ASSERT_THROW(log(A), IsNotZeroError);
  else
    ASSERT_NO_THROW(log(A));
}


TEST(BGptSU3Test, log){
  ptSU3 A(bgf::unit<bgf::AbelianBgf>());
  A.randomize();
  A.bgf().set_to_one();
  // this assumes A has unit background field
  ptsu3 B = log(A);
  // one loop
  ASSERT_TRUE( SU3Cmp(A[0], B[0])() );
  // two loop
  SU3 e = A[1] - A[0]*A[0]*0.5;
  ASSERT_TRUE( SU3Cmp(e, B[1])() );
  // three loop
  e = A[2] - (A[0]*A[1] + A[1]*A[0])*0.5 
    + A[0]*A[0]*A[0]/3;
  ASSERT_TRUE( SU3Cmp(e, B[2])() );
}

// Now, we want to check if log(exp(A)) = A holds...
TEST(BGptSU3Test, logexp){
  ptt::PtMatrix<MYORD> A = ptt::get_random_pt_matrix<MYORD>();
  ptt::PtMatrix<MYORD> B;
  B = log(exp<bgf::AbelianBgf,MYORD>(A));
  ptSU3 C{bgf::unit<bgf::AbelianBgf>()};
  ptSU3 D;
  C.randomize();
  C.bgf().set_to_one();
  D = exp<bgf::AbelianBgf,MYORD>(log(C));
  // higher orders
  // a quick and dirty estimate gives that the round-off error should
  // be around 100 * epsilon * n at order n ...
  double eps = std::numeric_limits<double>::epsilon() * 100;
  for (int r = 0; r < 6; ++r){
    ASSERT_TRUE( SU3Cmp(A[r], B[r])(eps * (r+1)) );
    ASSERT_TRUE( SU3Cmp(C[r], D[r])(eps * (r+1)) );
  }
}

// we need a main function here because we have to initialize the
// background field class...

int main(int argc, char **argv) {
  bgf::get_abelian_bgf(0,0,T, L);
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
