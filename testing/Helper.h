#ifndef _HELPER_H_
#define _HELPER_H_
#include "MyMath.h"
#include <limits>
#include <iostream>

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
  SU3Cmp(const SU3 &A, const SU3 &B) : a(A), b(B) {}
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
  (const double& eps = std::numeric_limits<double>::epsilon()*3)
    const {
    bool result = true;
    for (Cplx const * aptr = a.whr, * bptr = b.whr;
	 aptr != a.whr + 9 && result; ++aptr, ++bptr)
      if (fabs(aptr->re - bptr->re) > eps ||
	  fabs(aptr->im - bptr->im) > eps){
        std::cout << "at i = " << aptr - a.whr << std::endl;
        std::cout << aptr->re << ", " << aptr->im 
                  << std::endl
                  << bptr->re << ", " << bptr->im 
                  << std:: endl;
	result = false;
      }
    return result;
  }
  const SU3 &a, &b;
};
#endif
