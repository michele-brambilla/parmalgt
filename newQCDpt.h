#ifndef _MY_QCD_H_
#define _MY_QCD_H_

#include "MyQCD.h"
#include <algorithm>
#include <numeric>
#include <MyRand.h>
#include <iostream>
#include <Types.h>
#include <Background.h>


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
///  pt_q represents a series g q^(1) + g^2 q^(2) + ...
///
///  I decided this should be a separate class (from the U field) for
///  the following reason: If we multiply two q fields, the minimum
///  order that has to be considered decreases by one. Thus, one can
///  save time on the multiplications.
///  This class is (for now) only needed to be able to implement the
///  exp and log functions in a neat way...
///
///  \author Dirk Hesse <herr.dirk.hesse@gmail.com>
///  \date Tue Jan 24 16:55:10 2012

template <int AL_ORD, int PT_ORD> class pt_q {
public:
  // array type using template typedef
  typedef typename array_t<SU3, AL_ORD>::Type su3_array_t;

  pt_q() : min_ord_(0) { }

  SU3& operator[](const int& i) { return ptU_[i]; }
  const SU3& operator[](const int& i) const { return ptU_[i]; }

  su3_array_t& ptq() { return ptU_; };
  const su3_array_t& ptq() const { return ptU_; };

  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  ///
  ///  Multiplication-Assign.
  ///
  ///  Multiply each order by a given factor.
  ///
  ///  \author Dirk Hesse <herr.dirk.hesse@gmail.com>
  ///  \date Wed Jan 11 18:41:26 2012

  template <class C> pt_q& operator*=(const C &z) {
    for(int i = min_ord_; i < PT_ORD; i++)  ptU_[i] *= z;
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
  template <class C> pt_q& operator/=(const C &z) {
    for(int i = min_ord_; i < PT_ORD; i++)  ptU_[i] /= z;
    return *this;
  };


  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  ///
  ///  Addition-(Subtraction-)Assign for pt_q.
  ///
  ///  Here, we have to perform the operation order by order.
  ///
  ///  \author Dirk Hesse <herr.dirk.hesse@gmail.com>
  ///  \date Wed Jan 11 18:44:39 2012
  pt_q& operator+=(const pt_q &z) {
    for(int i = std::min(min_ord_, z.min_ord_); 
        i < PT_ORD; i++)  ptU_[i] += z[i];
    return *this;
  };
  pt_q& operator-=(const pt_q &z) {
    for(int i = std::min(min_ord_, z.min_ord_);
        i < PT_ORD; i++)  ptU_[i] -= z[i];
    return *this;
  };
  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  ///
  ///  Multipy-Assign.
  ///
  ///  This version is specialized for the pt_q type.
  ///
  ///  \author Dirk Hesse <herr.dirk.hesse@gmail.com>
  ///  \date Wed Jan 11 18:45:50 2012
  pt_q& operator*=(const pt_q& A) {
    su3_array_t tmp;
    // Handle products that involve pert. orders > 1
    for (int i = A.min_ord_; i < PT_ORD - 1; ++i)
      for (int j = min_ord_; j < PT_ORD  - 1 - i; ++j)
	tmp[i + j + 1] += ptU_[j] * A[i];
    // copy
    std::swap(tmp, ptU_);
    min_ord_ += A.min_ord_ + 1;
    return *this;
  };

  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  ///
  ///  Randiomize the matrices for each order, not the backgorund
  ///
  ///  \author Dirk Hesse <herr.dirk.hesse@gmail.com>
  ///  \date Wed Jan 11 18:46:17 2012

  void randomize() {
    static MyRand r(1235431);
    for (int i = 0; i < PT_ORD; ++i)
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

  template<class C> pt_q operator*(const C& z) const {
    pt_q result(*this);
    return result *= z;
  };
  template <class C> pt_q operator/(const C& z) const{
    return *this * (1./z);
  };
  template <class C> pt_q operator+(const C& z) const{
    pt_q result(*this);
    return result += z;
  };
  template <class C> pt_q operator-(const C& z) const{
    pt_q result(*this);
    return result -= z;
  };
private:
  su3_array_t ptU_;
  // minimum order to be taken into account on multiplication
  int min_ord_;
};


//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
///
///  BGptSU3 represents a series  V + g U^(1) + g^2 U^(2) + ...
///
///  \author Dirk Hesse <herr.dirk.hesse@gmail.com>
///  \date Tue Jan 24 16:55:46 2012

template <class B, int AL_ORD, int PT_ORD>
class BGptSU3{
 public:
  // array type using template typedef
  typedef typename array_t<SU3, AL_ORD>::Type su3_array_t;

  explicit BGptSU3(const B& bgf) : bgf_(bgf) { }
  // Default constructor, needed such that BGptSU3 may be used in a
  // std::array type.
  BGptSU3() : bgf_() { }
  SU3& operator[](const int& i) { return ptU_[i]; }
  const SU3& operator[](const int& i) const { return ptU_[i]; }
  B& bgf() { return bgf_; }
  const B& bgf() const { return bgf_; }
  
  su3_array_t& ptU() { return ptU_; };
  const su3_array_t& ptU() const { return ptU_; };

  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  ///
  ///  Multiplication-Assign.
  ///
  ///  Multiply each order by a given factor.
  ///
  ///  \author Dirk Hesse <herr.dirk.hesse@gmail.com>
  ///  \date Wed Jan 11 18:41:26 2012

  template <class C> BGptSU3& operator*=(const C &z) {
    for(int i = 0; i < PT_ORD; i++)  ptU_[i] *= z;
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
  template <class C> BGptSU3& operator/=(const C &z) {
    for(int i = 0; i < PT_ORD; i++)  ptU_[i] /= z;
    bgf_ /= z;
    return *this;
  };

  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  ///
  ///  Addition-Assign
  ///
  ///  Adding a scalar will only affect the lowest order in the
  ///  expansion.
  ///
  ///  \author Dirk Hesse <herr.dirk.hesse@gmail.com>
  ///  \date Wed Jan 11 18:43:19 2012
  template <class C> BGptSU3& operator+=(const C &z) {
    bgf_ += z;
    return *this;
  };
  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  ///
  ///  Subtract-Assign.
  ///
  ///  c.f. operator -=
  ///
  ///  \author Dirk Hesse <herr.dirk.hesse@gmail.com>
  ///  \date Wed Jan 11 18:44:12 2012
  template <class C> BGptSU3& operator-=(const C &z) {
    bgf_ -= z;
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
  BGptSU3& operator+=(const BGptSU3 &z) {
    for(int i = 0; i < PT_ORD; i++)  ptU_[i] += z[i];
    bgf_ += z.bgf();
    return *this;
  };
  BGptSU3& operator-=(const BGptSU3 &z) {
    for(int i = 0; i < PT_ORD; i++)  ptU_[i] -= z[i];
    bgf_ -= z.bgf();
    return *this;
  };
  BGptSU3& operator+=(const pt_q<AL_ORD, PT_ORD> &z) {
    for(int i = 0; i < PT_ORD; i++)  ptU_[i] += z[i];
    return *this;
  };
  BGptSU3& operator-=(const pt_q<AL_ORD, PT_ORD> &z) {
    for(int i = 0; i < PT_ORD; i++)  ptU_[i] -= z[i];
    return *this;
  };
  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  ///
  ///  Multipy-Assign.
  ///
  ///  This version is specialized for the BGptSU3 type.
  ///
  ///  \author Dirk Hesse <herr.dirk.hesse@gmail.com>
  ///  \date Wed Jan 11 18:45:50 2012
  BGptSU3& operator*=(const BGptSU3& A) {
    //SU3 tmp[allocORD];
    su3_array_t tmp;
    // Handle products that involve pert. orders > 1
    for (int i = 0; i < PT_ORD - 1; ++i)
      for (int j = 0; j < PT_ORD  - 1 - i; ++j)
	tmp[i + j + 1] += ptU_[j ] *A[i];
    // Include first order contributions
    for (int i = 0; i < PT_ORD; ++i)
      tmp[i] += bgf_.ApplyFromLeft(A[i]) +
	A.bgf_.ApplyFromRight(ptU_[i]);
    // Zeroth order.
    bgf_ *= A.bgf_;
    // copy
    std::swap(tmp, ptU_);
    return *this;
  };

  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  ///
  ///  Randiomize the matrices for each order, not the backgorund
  ///
  ///  \author Dirk Hesse <herr.dirk.hesse@gmail.com>
  ///  \date Wed Jan 11 18:46:17 2012

  void randomize() {
    static MyRand r(1235431);
    for (int i = 0; i < PT_ORD; ++i)
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

  template<class C> BGptSU3 operator*(const C& z) const {
    BGptSU3 result(*this);
    return result *= z;
  };
  template <class C> BGptSU3 operator/(const C& z) const{
    return *this * (1./z);
  };
  template <class C> BGptSU3 operator+(const C& z) const{
    BGptSU3 result(*this);
    return result += z;
  };
  template <class C> BGptSU3 operator-(const C& z) const{
    BGptSU3 result(*this);
    return result -= z;
  };

private:
  B bgf_;
  su3_array_t ptU_;
};

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
///
///  Dagger (hermitian conjugate) for ptSU3
///
///  \author Dirk Hesse <herr.dirk.hesse@gmail.com>
///  \date Thu Jan 26 11:47:34 2012

template <class B, int AL_ORD, int PT_ORD>
inline BGptSU3<B, AL_ORD, PT_ORD> 
dag( const BGptSU3<B, AL_ORD, PT_ORD>& U ){
  // The constructor call copies everyting
  BGptSU3<B, AL_ORD, PT_ORD> res(U);
  for (int i = 0; i < PT_ORD; ++i){
    // transpose
    std::swap(res[i][1], res[i][3]);
    std::swap(res[i][2], res[i][6]);
    std::swap(res[i][5], res[i][7]);
    // conjugate
    for (int j = 0; j < 9; ++j) res[i][j].im = -res[i][j].im;
  }
  return res;
};

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
///
///  Exponential
///
///  \author Dirk Hesse <herr.dirk.hesse@gmail.com>
///  \date Tue Jan 24 16:56:39 2012

template <class B, int AL_ORD, int PT_ORD>
inline BGptSU3<B, AL_ORD, PT_ORD> 
exp( const BGptSU3<B, AL_ORD, PT_ORD>& U ){
  // Make a copy of U (since exp U = 1 + U + ...)
  BGptSU3<B, AL_ORD, PT_ORD> result(U);
  // The pt orders are stored here:
  pt_q<AL_ORD, PT_ORD> q;
  // set q = U - V
  std::copy(result.ptU().begin(), result.ptU().end(),
            q.ptq().begin());
  pt_q<AL_ORD, PT_ORD> tmp(q); // store q^i here
  // calculate \tilde U = exp{q} up to order PT_ORD
  for (int i = 2; i <= PT_ORD; ++i){
    tmp *= q; // construct
    tmp /= i; // q^i / ( i! )
    result += tmp; // sum it up!
  }
  // convert to U^(i) = \tilde U^(i) V 
  for (int i = 0; i < PT_ORD; ++i)
    result[i] = result.bgf().ApplyFromRight(result[i]);
  return result;
}; 

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
///
///  Natural Logarithm  -- Using the Mercator series
///
///  \author Dirk Hesse <herr.dirk.hesse@gmail.com>
///  \date Tue Jan 24 16:56:55 2012

template <class B, int AL_ORD, int PT_ORD>
inline BGptSU3<B, AL_ORD, PT_ORD> 
log( const BGptSU3<B, AL_ORD, PT_ORD>& U ){
  // The constructor call copies the background field
  BGptSU3<B, AL_ORD, PT_ORD> result;
  result.bgf() = U.bgf();
  // The pt orders are stored here:
  pt_q<AL_ORD, PT_ORD> Util_m_one; // \tilde U - 1
  // convert to \tilde U^(i) = U^(i) V^{-1}
  B Vinv = U.bgf().inverse();
  for (int i = 0; i < PT_ORD; ++i){
    Util_m_one[i] = Vinv.ApplyFromRight(U[i]);
    result[i] = Util_m_one[i]; // now result = \tilde U - 1
  }
  pt_q<AL_ORD, PT_ORD> tmp(Util_m_one);
  // calculate q = log(\tilde U) up to order PT_ORD
  int sign = 1;
  for (int i = 2; i <= PT_ORD; ++i){
    sign *= -1;
    tmp *= Util_m_one; // construct Util^i
    result += tmp/i*sign; // sum it up!
  }
  return result;
};


template <class BGF,  int AL_ORD, int PT_ORD, int DIM>
class BGptGluon {
public:
  typedef BGptSU3<BGF, AL_ORD, PT_ORD> pt_su3_t;
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
    for (iterator i = begin(); i != end(); ++i) *i *= other;
    return *this;
  }

  template <class C>
  self_t& operator/=(const C& other){
    for (iterator i = begin(); i != end(); ++i) *i /= other;
    //for ( auto i : *this ) *i /= other;
    return *this;
  }

  pt_su3_t operator*(const self_t& other) const{
    return std::inner_product(begin(), end(),
                              other.begin(), pt_su3_t());
  }
  
private:
  array_t U_;
};

template  <int AL_ORD, int PT_ORD>
class BGptCVector {
public:
  typedef typename array_t<CVector, AL_ORD + 1>::Type vec_t;
  typedef BGptCVector self_t;
  typedef typename vec_t::iterator iterator;
  typedef typename vec_t::const_iterator const_iterator;
  // template tyedef for multiplicaion wit BGptSU3
  template <class BGF> struct pt_su3_t {
    typedef BGptSU3<BGF, AL_ORD, PT_ORD> Type;
  };

  /// Iterators
  iterator begin() { return v_.begin(); }
  iterator end() { return v_.end(); }
  const_iterator begin() const { return v_.begin(); }
  const_iterator end() const { return v_.end(); }
  /// Iterators
  iterator pt_end() { return v_.begin() + PT_ORD + 1; }
  const_iterator pt_end() const { return v_.begin() + PT_ORD + 1; }

  /// Access operators
  CVector& operator[](const int& i){ return v_[i]; }
  const CVector& operator[](const int& i) const { return v_[i]; }

  /// Arithmetic
  self_t& operator+=(const self_t& other){
    const_iterator j = other.begin();
    for (iterator i = begin(); 
         i != pt_end(); ++i, ++j) *i += *j;
    return *this;
  }
  self_t& operator-=(const self_t& other){
    const_iterator j = other.begin();
    for (iterator i = begin(); 
         i != pt_end(); ++i, ++j) *i += *j;
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
    for(iterator i = begin(); i != pt_end(); ++i)
      *i *= alpha;
    return *this;
  }
  template <class C>
  self_t& operator/=(const C& alpha) {
    for(iterator i = begin(); i != pt_end(); ++i)
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
  self_t operator*(const BGptSU3<BGF, AL_ORD, PT_ORD>& other) const{
    self_t result;
    // Handle products that involve pert. orders > 1
    // note that other[i] is of perturbative order i+1!
    for (int i = 1; i < PT_ORD + 1; ++i)
      for (int j = 0; j < i; ++j)
        result[i] += v_[j] * other[i - j - 1];
    
    // Handle the products involving other.bgf
    for (int i = 0; i < PT_ORD + 1; ++i)
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

template  <int AL_ORD, int PT_ORD, int DIM>
class BGptSpinColor {
public:
  typedef BGptCVector<AL_ORD, PT_ORD> pt_vec_t;
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
  self_t operator*(const BGptGluon<BGF, AL_ORD, PT_ORD, DIM>& U) 
    const {
    self_t result;
    // outer iteration 0 ... DIM
    const_iterator this_mu = begin();
    iterator psi_mu = result.begin();
    typename BGptGluon<BGF, AL_ORD, PT_ORD, DIM>\
      ::const_iterator U_mu = U.begin();
    for(; this_mu != end(); ++this_mu, ++psi_mu, ++U_mu)
      *psi_mu = *this_mu* *U_mu;
    return result;
  }
private:
  vec_t psi_;
  
};

#endif
