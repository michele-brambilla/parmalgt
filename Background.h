#ifndef _BACKGROUND_H_
#define _BACKGROUND_H_

#include "MyMath.h"
#include <vector>
#include <algorithm>
#include <functional>
#include <math.h>
#include <Types.h>

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
/// Background field classes.
///
/// Classes to implement trivial and abelian background fields.

namespace bgf {

  template <class C> inline SU3 operator*(const C& x, const SU3& y){
    return x.ApplyFromLeft(y);
  }
  template <class C> inline SU3 operator*( const SU3& y, const C& x){
    return x.ApplyFromRight(y);
  }
  
  /// Base class to define interface.
  class BgfBase {
  public:
    /// Left multiplication.
    virtual SU3 ApplyFromLeft ( const SU3 & ) const = 0;
    virtual CVector ApplyFromLeft ( const CVector& ) const = 0;
    /// Right multiplication.
    virtual SU3 ApplyFromRight ( const SU3 & ) const = 0;
    virtual CVector ApplyFromRight ( const CVector & ) const = 0;
  };

  /// Trivial (unit) background field.
  class TrivialBgf : public BgfBase {
  public:
    virtual SU3 ApplyFromLeft ( const SU3 & U) const {
      return U;
    }
    virtual SU3 ApplyFromRight ( const SU3 & U) const {
      return U;
    }
    virtual CVector ApplyFromLeft ( const CVector & U) const {
      return U;
    }
    virtual CVector ApplyFromRight ( const CVector & U) const {
      return U;
    }
    template <class C> TrivialBgf & operator*= (const C&) {
      return *this;
    }
    template <class C> TrivialBgf operator* (const C&) const {
      return *this;
    }
  };

  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  ///
  ///  Abelian \f$SU(3)\f$ background field.
  ///
  ///  \author Dirk Hesse <herr.dirk.hesse@gmail.com>
  ///  \date Tue Sep 27 11:08:30 2011
  class AbelianBgf : public BgfBase {
  public:

    //////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////
    ///
    ///  Vector constructor.
    ///
    ///  Construct from a three_vec_t as the diagonal elements.
    ///
    ///  \param v Vector we want to use for initialzing.
    ///  \author Dirk Hesse <herr.dirk.hesse@gmail.com>
    ///  \date Tue Sep 27 11:07:14 2011
    explicit AbelianBgf(const three_vec_t &v) : v_(v){ }
    AbelianBgf() : v_(){ }
    Cplx & operator[](const short& s){
      return v_[s];
    }
    //////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////
    ///
    ///  Left multiplication with an \f$SU(3)\f$ matrix.
    ///
    ///  This returns the product \f$ V U \f$, where \f$ V \f$ is the
    ///  background field's value represented by the instance.
    ///
    ///  \param U The \$f SU3 \$f matrix the instance is to be applied
    ///  to. 
    ///  \author Dirk Hesse <herr.dirk.hesse@gmail.com>
    ///  \date Tue Sep 27 11:03:50 2011
    virtual SU3 ApplyFromLeft ( const SU3 & U) const {
      SU3 result;
      for (int i = 0; i < 3; ++i)
	for (int j = 0; j < 3; ++j)
	  result.whr[3*i + j] = v_[i] * U.whr[3*i + j];
      return result;
    }
    virtual CVector ApplyFromLeft ( const CVector & v) const {
      CVector result;
      for (int i = 0; i < 3; ++i)
        result.whr[i] = v_[i] * v.whr[i];
      return result;
    }

    //////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////
    ///
    ///  Right multiplication with an \f$SU(3)\f$ matrix.
    ///
    ///  This returns the product \f$ U V \f$, where \f$ V \f$ is the
    ///  background field's value represented by the instance.
    ///
    ///  \param U The \$f SU3 \$f matrix the instance is to be applied
    ///  to. 
    ///  \author Dirk Hesse <herr.dirk.hesse@gmail.com>
    ///  \date Tue Sep 27 11:03:50 2011
    virtual SU3 ApplyFromRight ( const SU3 & U) const {
      SU3 result;
      for (int i = 0; i < 3; ++i)
	for (int j = 0; j < 3; ++j)
	  result.whr[3*i + j] = v_[j] * U.whr[3*i + j];
      return result;
    }
    virtual CVector ApplyFromRight ( const CVector & v) const {
      CVector result;
      for (int i = 0; i < 3; ++i)
        result.whr[i] = v_[i] * v.whr[i];
      return result;
    }
    bool operator==(const AbelianBgf& other) const{
      return v_ == other.v_;
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
		       std::bind1st( std::multiplies<Cplx>(), alpha));
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
    
    //////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////
    ///
    ///  Addition and subtraction of a sclar.
    ///
    ///  \author Dirk Hesse <herr.dirk.hesse@gmail.com>
    ///  \date Wed Jan 11 18:03:06 2012
    template<class C>
    AbelianBgf& operator+= (const C& alpha ){
      std::transform ( v_.begin(), v_.end(), v_.begin(),
		       std::bind1st( std::plus<C>(), alpha));
      return *this;
    }
    AbelianBgf& operator+= (const AbelianBgf& alpha ){
      std::transform (v_.begin(), v_.end(), alpha.v_.begin(),
                      v_.begin(), std::plus<Cplx>());
      return *this;
    }
    template<class C>
    AbelianBgf& operator-= (const C& alpha ){
      return *this += (-alpha);
    }
    template<class C>
    AbelianBgf operator+ (const C& alpha ) const {
      AbelianBgf result(*this);
      return result += alpha;
    }
    template<class C>
    AbelianBgf operator- (const C& alpha ) const {
      AbelianBgf result(*this);
      return result -= alpha;
    }
  private:
    //////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////
    ///
    ///  Diagonal elements of V.
    ///
    ///  \author Dirk Hesse <herr.dirk.hesse@gmail.com>
    ///  \date Tue Sep 27 10:50:29 2011
    three_vec_t v_;
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
 
  class AbelianBgfFactory {
  public:
    //////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////
    ///
    ///  Initialization of Vt_.
    ///
    ///  Here, we construct \f$ V(t) \f$ a la Peter Weisz, viz.
    ///  \f[ V(t) = \exp i ( {\cal E} x_0 - iC) \f],
    ///  where
    ///  \f[ {\cal E} = -i (C' - C) / T, \f]
    ///  and $C',C$ are the gluon boundary fields.
    ///
    ///  \param T Temporal lattice extend.
    ///  \param L Spatial lattice extend.
    ///  \param eta Background field parameter. C.f. P. Weisz. Usually
    ///  0 is used.
    ///  \param nu Background field parameter. C.f. P. Weisz. Usually
    ///  0 is used.
    ///
    ///  \author Dirk Hesse <herr.dirk.hesse@gmail.com>
    ///  \date Tue Sep 27 10:46:16 2011
    AbelianBgfFactory(const int &T, const int& L,
                      const double& eta = 0.0, const double& nu = 0.0) 
      : Vt_(T){
      if (T < 0)
        throw std::exception();
      double pi = std::atan(1.)*4.;
      double gamma = 1./L/T * (eta + pi/3.);
      double eps[3] = {-2.*gamma, gamma, gamma};
      double iC[3] = {-( eta - pi/3) / L,
		      -eta * (-0.5 + nu)  / L,
		      -( eta * (0.5 + nu) + pi/3.) / L};
      for (int t = 0; t < T; ++t)
	for (int k = 0; k < 3; ++k)
	  Vt_[t][k] = exp(Cplx(0,eps[k] * t - iC[k]));
    };
    const three_vec_t& get(int t){ return Vt_.at(t); }
  private:
    std::vector <three_vec_t> Vt_;
  };

  inline AbelianBgf get_abelian_bgf(const int& t,
                                    const int& mu,
                                    const int& T = -1,
                                    const int& L = 0,
                                    const double& eta = 0.0,
                                    const double& nu = 0.0){
    static AbelianBgfFactory factory(T, L, eta, nu);
    static three_vec_t one = {1,1,1};
    if (mu)
      return AbelianBgf(one);
    else
      return AbelianBgf(factory.get(t));
  };
  
  inline AbelianBgf unit(const Cplx& alpha = Cplx(1.0,0.)){
    three_vec_t alpha_v = {alpha, alpha, alpha};
    return AbelianBgf(alpha_v);
  };
}

#endif
