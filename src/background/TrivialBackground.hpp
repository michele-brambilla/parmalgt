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
  template<class Matrix_t = SU3, class Vector_t = CVector>
  class TrivialBgf : public BgfBase {
  public:
    Matrix_t ApplyFromLeft ( const Matrix_t & U) const override {
      return U;
    }
    Matrix_t ApplyFromRight ( const Matrix_t & U) const override {
      return U;
    }
    Matrix_t Add (const Matrix_t & U) const override {
      Matrix_t res(U);
      for ( int i = 0; i < Matrix_t::size ;++i) res(i,i) += 1;
      return res;
    }
    Vector_t ApplyFromLeft ( const Vector_t & U) const override {
      return U;
    }
    Vector_t ApplyFromRight ( const Vector_t & U) const override {
      return U;
    }
    template <class C> TrivialBgf & operator*= (const C&) {
      return *this;
    }
    template <class C> TrivialBgf operator* (const C&) const {
      return *this;
    }
    TrivialBgf inverse() const { return *this; }
    void set_to_one() override { }
    void set_to_zero() override { throw std::exception(); } // not possible ...
    TrivialBgf dag() const { return *this; }
    Cplx Tr() const override { return 3; }
    double Norm() const override { return 1; }
    void Trless() override {}
    void reH() override {}
  };


}

#endif // _TRIVIAL_BACKGROUND_H_
