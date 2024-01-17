#ifndef _TRIVIAL_BACKGROUND_H_
#define _TRIVIAL_BACKGROUND_H_

#include "BackgroundInterface.hpp"

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

  /// Trivial (unit) background field.
  class TrivialBgf /*: public BgfBase*/ {
  public:
    using value_type = double;
    static const size_t storage_size = 0;

    template<class T>
    T ApplyFromLeft ( const T & U) const {
      return U;
    }
    template<class T>
    T ApplyFromRight ( const T& U) const {
      return U;
    }
    SU3 Add (const SU3 & U) const {
      SU3 res(U);
      for ( int i = 0; i < SU3::size ;++i) res(i,i) += 1;
      return res;
    }
    template <class C> TrivialBgf & operator*= (const C&) {
      return *this;
    }
    template <class C> TrivialBgf operator* (const C&) const {
      return *this;
    }
    TrivialBgf inverse() const { return *this; }
    void set_to_one() { }
    void set_to_zero() { throw std::exception(); } // not possible ...
    TrivialBgf dag() const { return *this; }
    Cplx Tr() const { return 3; }
    double Norm() const { return 1; }
    void Trless() {}
    void reH() {}
  };


}

#endif // _TRIVIAL_BACKGROUND_H_
