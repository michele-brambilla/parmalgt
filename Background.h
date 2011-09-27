
#ifndef _BACKGROUND_H_
#define _BACKGROUND_H_

#include "MyMath.h"
#include <vector>
#include <algorithm>
#include <functional>
#include <math.h>

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
    /// Right multiplication.
    virtual SU3 ApplyFromRight (const SU3 & ) const = 0;
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
    template <class C> TrivialBgf & operator*= (const C&) {
      return *this;
    }
    template <class C> TrivialBgf operator* (const C&) const {
      return *this;
    }
  };

  class InitializedWithWrongSize : public std::exception { };

  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  ///
  ///  Abelian \f$SU3\f$ background field.
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
    ///  Construct from a vector of size 3 as the diagonal elements.
    ///
    ///  \param v Vector we want to use for initialzing.
    ///  \throws InitializedWithWrongSize Thrown if the length of v is
    ///  not equal to three.
    ///  \author Dirk Hesse <herr.dirk.hesse@gmail.com>
    ///  \date Tue Sep 27 11:07:14 2011
    explicit AbelianBgf(const std::vector<Cplx> &v) throw(InitializedWithWrongSize): v_(v){
      if (v_.size() != 3)
	throw InitializedWithWrongSize();
    }
    explicit AbelianBgf(const int &t) : v_(Vt_[t]) { };
    AbelianBgf() : v_(3){ }
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
    static void init(const int &T, const int& L,
		     const double& eta = 0.0, const double& nu = 0.0) {
      double pi = std::atan(1.)*4.;
      double gamma = 1./L/T * (eta + pi/3.);
      double eps[3] = {-2.*gamma, gamma, gamma};
      double iC[3] = {-( eta - pi/3) / L,
		      -eta * (-0.5 + nu)  / L,
		      -( eta * (0.5 + nu) + pi/3.) / L};
      Vt_ = std::vector< std::vector<Cplx> >(T, std::vector<Cplx>(3));
      for (int t = 1; t < T; ++t)
	for (int k = 0; k < 3; ++k)
	  Vt_[t][k] = exp(Cplx(0,eps[k] * t - iC[k]));
      for (int k = 0; k < 3; ++k)
	Vt_[0][k] = Cplx(1.,0.);
    }
  private:
    //////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////
    ///
    ///  Diagonal elements of V.
    ///
    ///  \author Dirk Hesse <herr.dirk.hesse@gmail.com>
    ///  \date Tue Sep 27 10:50:29 2011
    std::vector<Cplx> v_;
    //////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////
    ///
    ///  Background field. The Background field in the SF reads
    ///
    ///  \f[ V_\mu(x) = V_\mu(x_0), \quad V_0(x_0) = 1,\quad V_k(x_0)
    ///  = V(x_0), \f]
    ///
    ///  where each \f$ V(x_0) \f$ is a diagonal matrix. Vk_ stores
    ///  the diagonal elments of \f$ V(x_0) \f$. It must be
    ///  initialized with the init method before it can be used.
    ///
    ///  \author Dirk Hesse <herr.dirk.hesse@gmail.com>
    ///  \date Tue Sep 27 10:46:31 2011
    static std::vector <std::vector <Cplx> > Vt_;
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
