#ifndef _BACKGROUND_HELPER_HPP_
#define _BACKGROUND_HELPER_HPP_

#include <vector>
#include <algorithm>
#include <numeric>
#include <functional>
#include <math.h>
#include <iostream>

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
/// Background field classes.
///
/// Classes to implement trivial and abelian background fields.

namespace bgf {

  template <class Matrix_t, class C> inline Matrix_t operator*(const C& x, const Matrix_t& y){
    return x.ApplyFromLeft(y);
  }
  template <class Matrix_t, class C> inline Matrix_t operator*(const Matrix_t& y, const C& x){
    return x.ApplyFromRight(y);
  }

}

#endif // _BACKGROUND_HELPER_HPP_