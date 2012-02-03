#ifndef _MY_QCD_H_
#define _MY_QCD_H_

#include "MyQCD.h"
#include <algorithm>
#include <numeric>
#include <MyRand.h>
#include <iostream>
#include <Types.h>
#include <Background.h>
#include <PtTypes.hpp>
#include <testing/Helper.h>


// debug flag

const bool do_debug = true;

// is thrown by log if called with a field with a non-trivial lowest
// order.

class IsNotOneError : public std::exception { };

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
///
///  This is a simple struct to do the check if a field has a trivial
///  lowest order.
///
///  The point is that
///
///  a) You should not use #defines
///  b) A macro would be a pain anyway, because it must know,
///     depending on the background field class we use, what exactly
///     "non trivial" means. Hence a template is the way to go.
///
///  \author Dirk Hesse <herr.dirk.hesse@gmail.com>
///  \date Fri Feb  3 12:52:30 2012


template <bool shall_i_debug, class B> struct IsOne {
  static const bool debug_on = false;
  static void check(const B&) { } // do not debug
};

// overload for shall_i_debug == true

template <class B> struct IsOne<true, B> {
  static const bool debug_on = true;
  static void check(const B& V) {
    static double eps = std::numeric_limits<double>::epsilon() * 10;
    if ( (V - bgf::unit<B>()).Norm() > eps)  throw IsNotOneError();
  }
};

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
///
///  This is just a stub of the original class (which is of course
///  saved in the git repo). Only the arhitmetic methods are
///  implemented, the rest will follow.
///
///  \author Dirk Hesse <herr.dirk.hesse@gmail.com>
///  \date Tue Oct 11 16:23:38 2011


//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
///
///  BGptSU3 represents a series  V + g U^(1) + g^2 U^(2) + ...
///
///  \author Dirk Hesse <herr.dirk.hesse@gmail.com>
///  \date Tue Jan 24 16:55:46 2012

template <class B, int ORD>
class BGptSU3{
public:
  // self type
  typedef BGptSU3 self_t;
  // underlying perturbative structure
  typedef typename ptt::PtMatrix<ORD> pt_matrix_t;
  // iterators
  typedef typename pt_matrix_t::iterator iterator;
  typedef typename pt_matrix_t::const_iterator const_iterator;
  // 3x3 matrix type
  typedef typename pt_matrix_t::SU3_t SU3_t;

  explicit BGptSU3(const B& bgf) : bgf_(bgf) { }
  BGptSU3(const B& bgf, const pt_matrix_t& ptU) : 
    bgf_(bgf), ptU_(ptU) { }

  // Default constructor, needed such that self_t may be used in a
  // std::array type.
  BGptSU3() : bgf_(bgf::unit<B>()) { }

  SU3_t& operator[](const int& i) { return ptU_[i]; }
  const SU3_t& operator[](const int& i) const { return ptU_[i]; }

  B& bgf() { return bgf_; }
  const B& bgf() const { return bgf_; }
  
  pt_matrix_t& ptU() { return ptU_; };
  const pt_matrix_t& ptU() const { return ptU_; };

  const_iterator begin() const { return ptU_.begin(); }
  const_iterator end() const { return ptU_.end(); }
  iterator begin() { return ptU_.begin(); }
  iterator end() { return ptU_.end(); }

  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  ///
  ///  Multiplication-Assign.
  ///
  ///  Multiply each order by a given factor.
  ///
  ///  \author Dirk Hesse <herr.dirk.hesse@gmail.com>
  ///  \date Wed Jan 11 18:41:26 2012

  template <class C> self_t& operator*=(const C &z) {
    ptU_.multiply_assign(z);
    bgf_ *= z;
    return *this;
  };
  
  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  ///
  ///  Division-Assign. 
  ///
  ///  See operator *=.
  ///
  ///  \author Dirk Hesse <herr.dirk.hesse@gmail.com>
  ///  \date Wed Jan 11 18:42:49 2012
  template <class C> self_t& operator/=(const C &z) {
    ptU_.multiply_divide(z);
    bgf_ /= z;
    return *this;
  };


  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  ///
  ///  Addition-(Subtraction-)Assign for BGptSU3.
  ///
  ///  Here, we have to perform the operation order by order.
  ///
  ///  \author Dirk Hesse <herr.dirk.hesse@gmail.com>
  ///  \date Wed Jan 11 18:44:39 2012
  self_t& operator+=(const self_t &z) {
    ptU_.add_assign(z.ptU_);
    bgf_ += z.bgf();
    return *this;
  };
  self_t& operator-=(const self_t &z) {
    ptU_.sub_assign(z.ptU_);
    bgf_ -= z.bgf();
    return *this;
  };

  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  ///
  ///  Multipy-Assign.
  ///
  ///  This version is specialized for the BGptSU3 type.
  ///
  ///  We implement this by writing
  ///    (V1 + U1)(V2 + U2) = V1*V2 + V1*U2 + U1*V2 + U2*U2
  ///
  ///  \author Dirk Hesse <herr.dirk.hesse@gmail.com>
  ///  \date Wed Jan 11 18:45:50 2012

  self_t& operator*=(const self_t& A) {
    // U1*U2
    pt_matrix_t tmp = ptt::multiply(ptU_, A.ptU_);
    // The V*U terms
    tmp.add_assign(ptt::multiply(bgf_, A.ptU_));
    tmp.add_assign(ptt::multiply(A.bgf_, ptU_));
    // V1*V2
    bgf_ *= A.bgf_;
    // keep in mind we assign-multiply here
    std::swap(tmp, ptU_);
    return *this;
  }


  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  ///
  ///  Randiomize the matrices for each order, not the backgorund
  ///
  ///  \author Dirk Hesse <herr.dirk.hesse@gmail.com>
  ///  \date Wed Jan 11 18:46:17 2012

  void randomize() {
    static MyRand r(1235431);
    for (int i = 0; i < ORD; ++i)
      for (int j = 0; j < 3; ++j)
        for (int k = 0; k < 3; ++k)
          ptU_[i](j, k) = Cplx(r.Rand(), r.Rand());
  };

  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  ///
  ///  Mulitplication, Division, Addition, Subtraction.
  ///
  ///  To implement those, we reduce, reuse, recycle the
  ///  assign-operators! Not that we don't have to care about whether
  ///  or not we deal with a scalar, because the assign operators
  ///  already take care of this.
  ///
  ///  \author Dirk Hesse <herr.dirk.hesse@gmail.com>
  ///  \date Wed Jan 11 18:46:45 2012

  template<class C> self_t operator*(const C& z) const {
    self_t result(*this);
    return result *= z;
  };
  template <class C> self_t operator/(const C& z) const{
    return *this * (1./z);
  };
  template <class C> self_t operator+(const C& z) const{
    self_t result(*this);
    return result += z;
  };
  template <class C> self_t operator-(const C& z) const{
    self_t result(*this);
    return result -= z;
  };
  void id(){
    std::fill(begin(),end(),SU3());
    bgf_.set_to_one();
  }
  void zero(){
    std::fill(begin(),end(),SU3());
    bgf_.set_to_zero();
  }
  /// Trace
  void Tr(Cplx *tt) const {
    for(int i = 0; i < PTORD; i++)
      tt[i+1] = ptU_[i].whr[0] + ptU_[i].whr[4] + ptU_[i].whr[8];
    tt[0] = bgf_.Tr();
  }
  /// Make traceless
  void Trless(){
    Cplx z;
    for(int i = 0; i < PTORD; i++){
      z = D3*(ptU_[i][0] + ptU_[i][4] + ptU_[i][8]);
      ptU_[i][0] -= z;
      ptU_[i][4] -= z;
      ptU_[i][8] -= z;
    }
    bgf_.set_to_zero();
  }
  self_t& reH(){
    //for(int i = 0; i < PTORD; i++) ptU_[i].reH();
    //bgf_.reH();
    Cplx tr;
    for(int i = 0; i < PTORD; i++){
      ptU_[i] -= dag(ptU_[i]);
      ptU_[i] *= .5;
      tr = ptU_[i].Tr()/3.;
      ptU_[i][0] -= tr;
      ptU_[i][4] -= tr;
      ptU_[i][8] -= tr;
    }
    //bgf_.set_to_zero();
    return *this;
  }

private:
  B bgf_;
  pt_matrix_t ptU_;
};

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
///
///  Dagger (hermitian conjugate) for ptSU3
///
///  \author Dirk Hesse <herr.dirk.hesse@gmail.com>
///  \date Thu Jan 26 11:47:34 2012

//template <class B, int AL_ORD, int PT_ORD>
//inline BGptSU3<B, AL_ORD, PT_ORD> 
//dag( const BGptSU3<B, AL_ORD, PT_ORD>& U ){
//  // The constructor call copies everyting
//  BGptSU3<B, AL_ORD, PT_ORD> res(U);
//  for (int i = 0; i < PT_ORD; ++i){
//    // transpose
//    std::swap(res[i][1], res[i][3]);
//    std::swap(res[i][2], res[i][6]);
//    std::swap(res[i][5], res[i][7]);
//    // conjugate
//    for (int j = 0; j < 9; ++j) res[i][j].im = -res[i][j].im;
//  }
//  res.bgf() = U.bgf().dag();
//  return res;
//};

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
///
///  Exponential
///
///  \author Dirk Hesse <herr.dirk.hesse@gmail.com>
///  \date Tue Jan 24 16:56:39 2012

//template <class B, int AL_ORD, int PT_ORD>
//BGptSU3<B, AL_ORD, PT_ORD> 
//exp( const BGptSU3<B, AL_ORD, PT_ORD>& );

template <class B, int ORD>
BGptSU3<B, ORD> exp(const ptt::PtMatrix<ORD>& q){
  BGptSU3<B, ORD> result(bgf::unit<B>(), q);
  // result = 1 + q 
  ptt::PtMatrix<ORD> tmp(q);
  for ( int i = 2; i < ORD; ++i){
    tmp = ptt::multiply(tmp, q);
    tmp.divide_assign((double)i);
    result.ptU().add_assign(tmp);
  }
  return result;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
///
///  Natural Logarithm  -- Using the Mercator series
///
///  FIXME: THIS does not work
///
///  \author Dirk Hesse <herr.dirk.hesse@gmail.com>
///  \date Tue Jan 24 16:56:55 2012

template <class B, int ORD>
inline ptt::PtMatrix<ORD> log(const BGptSU3<B, ORD>& U){
  IsOne<do_debug, B>::check(U.bgf());
  ptt::PtMatrix<ORD> result(U.ptU()), tmp(U.ptU());
  double sign = 1.;
  for (int i = 2; i <= ORD; ++i){
    sign *= -1;
    tmp = ptt::multiply(tmp, U.ptU());
    result.add_assign( ptt::multiply(tmp, (double)i/sign) ); 
    // FIXME stretch_assign would be good here!
  }
  return result;
};

//template <class B, int AL_ORD, int PT_ORD>
//inline BGptSU3<B, AL_ORD, PT_ORD> 
//log( const BGptSU3<B, AL_ORD, PT_ORD>& U ){
//  // The constructor call copies the background field
//  BGptSU3<B, AL_ORD, PT_ORD> result;
//  result.bgf() = U.bgf();
//  // The pt orders are stored here:
//  pt_q<AL_ORD, PT_ORD> Util_m_one; // \tilde U - 1
//  // convert to \tilde U^(i) = U^(i) V^{-1}
//  B Vinv = U.bgf().inverse();
//  for (int i = 0; i < PT_ORD; ++i){
//    Util_m_one[i] = Vinv.ApplyFromRight(U[i]);
//    result[i] = Util_m_one[i]; // now result = \tilde U - 1
//  }
//  pt_q<AL_ORD, PT_ORD> tmp(Util_m_one);
//  // calculate q = log(\tilde U) up to order PT_ORD
//  int sign = 1;
//  for (int i = 2; i <= PT_ORD; ++i){
//    sign *= -1;
//    tmp *= Util_m_one; // construct Util^i
//    result += tmp/i*sign; // sum it up!
//  }
//  return result;
//};


template <class BGF,  int ORD, int DIM>
class BGptGluon {
public:
  typedef BGptSU3<BGF, ORD> pt_su3_t;
  typedef typename array_t<pt_su3_t, DIM>::Type array_t;
  typedef BGptGluon self_t;
  typedef typename array_t::iterator iterator;
  typedef typename array_t::const_iterator const_iterator;
  // access
  pt_su3_t& operator[](const int& i){ return U_[i]; }
  const pt_su3_t& operator[](const int& i) const { return U_[i]; }
  // iterators
  iterator begin(){return U_.begin();}
  const_iterator begin() const {return U_.begin();}
  iterator end(){return U_.end();}
  const_iterator end() const {return U_.end();}

  template <class C>
  self_t& operator*=(const C& other){
#ifdef HAVE_CXX0X
    for (auto& e : U_) { e *= other; }
#else
    for (iterator i = begin(); i != end(); ++i) *i *= other;
#endif
    return *this;
  }

  template <class C>
  self_t& operator/=(const C& other){
 #ifdef HAVE_CXX0X
    for (auto& e : U_) { e /= other; }
#else
    for (iterator i = begin(); i != end(); ++i) *i /= other;
#endif
     return *this;
  }

  pt_su3_t operator*(const self_t& other) const{
    return std::inner_product(begin(), end(),
                              other.begin(), pt_su3_t());
  }
  
private:
  array_t U_;
};

template  <int ORD>
class BGptCVector {
public:
  typedef typename array_t<CVector, ORD + 1>::Type vec_t;
  typedef BGptCVector self_t;
  typedef typename vec_t::iterator iterator;
  typedef typename vec_t::const_iterator const_iterator;
  // template tyedef for multiplicaion wit BGptSU3
  template <class BGF> struct pt_su3_t {
    typedef BGptSU3<BGF, ORD> Type;
  };

  /// Iterators
  iterator begin() { return v_.begin(); }
  iterator end() { return v_.end(); }
  const_iterator begin() const { return v_.begin(); }
  const_iterator end() const { return v_.end(); }

  /// Access operators
  CVector& operator[](const int& i){ return v_[i]; }
  const CVector& operator[](const int& i) const { return v_[i]; }

  /// Arithmetic
  self_t& operator+=(const self_t& other){
    const_iterator j = other.begin();
    for (iterator i = begin(); 
         i != end(); ++i, ++j) *i += *j;
    return *this;
  }
  self_t& operator-=(const self_t& other){
    const_iterator j = other.begin();
    for (iterator i = begin(); 
         i != end(); ++i, ++j) *i += *j;
    return *this;
  }
  self_t operator+(const self_t& other) const{
    self_t result(*this);
    return result += other;
  }
  self_t operator-(const self_t& other) const{
    self_t result(*this);
    return result -= other;
  }
  template <class C>
  self_t& operator*=(const C& alpha) {
#ifdef HAVE_CXX0X
    for (auto& e : v_) { e *= alpha; }
#else
    for(iterator i = begin(); i != end(); ++i)
      *i *= alpha;
#endif
    return *this;
  }
  template <class C>
  self_t& operator/=(const C& alpha) {
#ifdef HAVE_CXX0X
    for (auto& e : v_) { e /= alpha; }
#else
    for(iterator i = begin(); i != end(); ++i)
      *i /= alpha;
#endif
    return *this;
  }
  template <class C>
  self_t operator*(const C& alpha) const {
    self_t result;
    return result *= alpha;
  }
  template <class C>
  self_t operator/(const C& alpha) const {
    self_t result;
    return result /= alpha;
  }  
  template <class BGF>
  self_t operator*(const BGptSU3<BGF, ORD>& other) const{
    self_t result;
    // Handle products that involve pert. orders > 1
    // note that other[i] is of perturbative order i+1!
    for (int i = 1; i < ORD + 1; ++i)
      for (int j = 0; j < i; ++j)
        result[i] += v_[j] * other[i - j - 1];
    
    // Handle the products involving other.bgf
    for (int i = 0; i < ORD + 1; ++i)
      result[i] += other.bgf().ApplyFromLeft(v_[i]);
    return result;
  }

  bool operator==(const self_t& other) const {
    for (const_iterator i = begin(), j = other.begin();
         i != end(); ++i, ++j)
      if (*i != *j)
        return false;
    return true;
  }
private:
  vec_t v_;
};

template  <int ORD, int DIM>
class BGptSpinColor {
public:
  typedef BGptCVector<ORD> pt_vec_t;
  typedef typename array_t<pt_vec_t, DIM>::Type vec_t;
  typedef typename vec_t::iterator iterator;
  typedef typename vec_t::const_iterator const_iterator;
  typedef BGptSpinColor self_t;
  
  // access
  pt_vec_t& operator[](const int& i){ return psi_[i]; }
  const pt_vec_t& operator[](const int& i) const { return psi_[i]; }
  
  // iterators
  iterator begin() { return psi_.begin(); }
  iterator end() { return psi_.end(); }
  const_iterator begin() const { return psi_.begin(); }
  const_iterator end() const { return psi_.end(); }

  // arithmetic with self
  self_t& operator+=(const self_t& other){
    const_iterator j = other.begin();
    for (iterator i = begin(); i != end(); ++i, ++j) *i += *j;
    return *this;
  }
  self_t& operator-=(const self_t& other){
    const_iterator j = other.begin();
    for (iterator i = begin(); i != end(); ++i, ++j) *i -= *j;
    return *this;
  }
  self_t operator+(const self_t& other) const{
    self_t result(*this);
    return result += other;
  }
  self_t operator-(const self_t& other) const{
    self_t result(*this);
    return result -= other;
  }
  // scalar arithmetic
    template <class C>
  self_t& operator*=(const C& alpha) {
    for(iterator i = begin(); i != end(); ++i)
      *i *= alpha;
    return *this;
  }
  template <class C>
  self_t& operator/=(const C& alpha) {
    for(iterator i = begin(); i != end(); ++i)
      *i /= alpha;
    return *this;
  }
  template <class C>
  self_t operator*(const C& alpha) const {
    self_t result;
    return result *= alpha;
  }
  template <class C>
  self_t operator/(const C& alpha) const {
    self_t result;
    return result /= alpha;
  }
  template <class BGF>
  self_t operator*(const BGptGluon<BGF, ORD, DIM>& U) 
    const {
    self_t result;
    // outer iteration 0 ... DIM
    const_iterator this_mu = begin();
    iterator psi_mu = result.begin();
    typename BGptGluon<BGF, ORD, DIM>
      ::const_iterator U_mu = U.begin();
    for(; this_mu != end(); ++this_mu, ++psi_mu, ++U_mu)
      *psi_mu = *this_mu* *U_mu;
    return result;
  }
  void uno_p_gmu(SpinColor&, int, int);
  void uno_m_gmu(SpinColor&, int, int);

private:
  vec_t psi_;
  
};

#endif
