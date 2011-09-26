#ifndef _BACKGROUND_H_
#define _BACKGROUND_H_

#include "MyMath.h" // SU3
#include <vector>
#include <algorithm>
#include <functional>

///
/// Background field classes

namespace bgf {

  template <class C> inline SU3 operator*(const C& x, const SU3& y){
    return x.ApplyFromLeft(y);
  }
  template <class C> inline SU3 operator*( const SU3& y, const C& x){
    return x.ApplyFromRight(y);
  }
  
  /// Base class to define interface
  class BgfBase {
  public:
    virtual SU3 ApplyFromLeft ( const SU3 & ) const = 0;
    virtual SU3 ApplyFromRight (const SU3 & ) const = 0;
  };

  class TrivialBgf : public BgfBase {
  public:
    virtual SU3 ApplyFromLeft ( const SU3 & U) const {
      return U;
    }
    virtual SU3 ApplyFromRight ( const SU3 & U) const {
      return U;
    }
    template <class C> TrivialBgf & operator*= (const C&) {
      return *this;
    }
    template <class C> TrivialBgf operator* (const C&) const {
      return *this;
    }
  };

  class InitializedWithWrongSize : public std::exception { };

  class AbelianBgf : public BgfBase {
  public:
    explicit AbelianBgf(const std::vector<Cplx> &v) : v_(v){
      if (v_.size() != 3)
	throw InitializedWithWrongSize();
    }
  AbelianBgf() : v_(3){ }
    Cplx & operator[](const short& s){
      return v_[s];
    }
    virtual SU3 ApplyFromLeft ( const SU3 & U) const {
      SU3 result;
      for (int i = 0; i < 3; ++i)
	for (int j = 0; j < 3; ++j)
	  result.whr[3*i + j] = v_[i] * U.whr[3*i + j];
      return result;
    }
    virtual SU3 ApplyFromRight ( const SU3 & U) const {
      SU3 result;
      for (int i = 0; i < 3; ++i)
	for (int j = 0; j < 3; ++j)
	  result.whr[3*i + j] = v_[j] * U.whr[3*i + j];
      return result;
    }
    //////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////
    ///
    ///  Generic *= operator template.
    ///
    ///  \author Dirk Hesse <herr.dirk.hesse@gmail.com>
    ///  \date Mon Sep 26 18:28:23 2011
    template<class C>
    AbelianBgf& operator*= ( const C& alpha ) {
      std::transform ( v_.begin(), v_.end(), v_.begin(),
		       std::bind1st( std::multiplies<C>(), alpha));
      return *this;
    }
    //////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////
    ///
    ///  Generic multiplication.
    ///
    ///  \author Dirk Hesse <herr.dirk.hesse@gmail.com>
    ///  \date Mon Sep 26 18:28:36 2011
    template<class C>
    AbelianBgf operator* (const C& alpha ) const {
      AbelianBgf result(*this);
      return result *= alpha;
    }
  private:
    std::vector<Cplx> v_;
  };

  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  ///
  ///  Specialization of *= and * operator for AbelianBgf type.
  ///
  ///  The Arithmetics are quite straight forward, altough the
  ///  AbelianBgf*AbelianBgf multiplication is a special case.
  ///
  ///  \author Dirk Hesse <herr.dirk.hesse@gmail.com>
  ///  \date Mon Sep 26 18:25:58 2011

  template<> inline AbelianBgf& AbelianBgf::operator*= 
  ( const AbelianBgf& other ){
    std::transform  (v_.begin(), v_.end(), other.v_.begin(),
		     v_.begin(), std::multiplies<Cplx>() );
    return *this;
  }
  
}

#endif
