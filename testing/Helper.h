#ifndef _HELPER_H_
#define _HELPER_H_
#include "MyMath.h"

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
#endif
