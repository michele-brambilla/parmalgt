#ifndef _BACKGROUND_H_
#define _BACKGROUND_H_

#include "MyMath.h" // SU3
#include <vector>

///
/// Background field classes

namespace bgf {

  template <class C> inline SU3 operator*(const C& x, const SU3& y){
    x.ApplyFromLeft(y);
  }
  template <class C> inline SU3 operator*( const SU3& y, const C& x){
    x.ApplyFromRight(y);
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
  };

  class AbelianBgf : public BgfBase {
  public:
    virtual SU3 ApplyFromLeft ( const SU3 & U) const {
      SU3 result;
      for (int i = 0; i < 3; ++i)
	for (int j = 0; j < 3 +=j)
	  result[3*i + j] = v_[j] * U[3*i + j]
    }
  private:
    std::vector<Cplx> v_;
  };
  
}

#endif
