#include "Background.h"
#include "gtest/gtest.h"

// We need a test for the AbelianBgf class...

// With a nice base class ...

MyRand r(23797);


//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
///
///  Class to compare SU3 Matrix types
///
///  \author Dirk Hesse <herr.dirk.hesse@gmail.com>
///  \date Wed Sep 21 19:12:10 2011
struct SU3Cmp {
  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  ///
  ///  Constructor.
  ///
  ///  \param A First of the matrices to be compared
  ///  \param B Second matrix to be compared 
  ///
  ///  \author Dirk Hesse <herr.dirk.hesse@gmail.com>
  ///  \date Wed Sep 21 19:12:29 2011
  SU3Cmp(SU3 A, SU3 B) : a(A), b(B) {}
  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  ///
  ///  Perform the Comparison
  ///
  ///  \param eps Accuracy to which we want to the matrices' entries
  ///  to be compared. Default value is three machine accuracies.
  ///
  ///  \author Dirk Hesse <herr.dirk.hesse@gmail.com>
  ///  \date Wed Sep 21 19:14:08 2011
  bool operator()
  (const double& eps = std::numeric_limits<double>::epsilon()*3){
    bool result = true;
    for (Cplx *aptr = a.whr, *bptr = b.whr;
	 aptr != a.whr + 9 && result;
	 ++aptr, ++bptr)
      if (fabs(aptr->re - bptr->re) > eps ||
	  fabs(aptr->im - bptr->im) > eps)
	result = false;
    return result;
  }
  SU3 &a, &b;
};

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
