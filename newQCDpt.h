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
    //ptU_ = tmp;
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
