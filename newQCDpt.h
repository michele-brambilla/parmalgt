#ifndef _MY_QCD_H_
#define _MY_QCD_H_

#include "MyQCD.h"
#include <algorithm>
#include <MyRand.h>
#include <iostream>

class ptSU3;
class ptCVector;
class ptGluon;
class ptSpinColor;

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
///
///  This is just a stub of the original class (which is of course
///  saved in the git repo). Only the arhitmetic methods are
///  implemented, the rest will follow.
///
///  \author Dirk Hesse <herr.dirk.hesse@gmail.com>
///  \date Tue Oct 11 16:23:38 2011

template <class B>
class BGptSU3{
 private:
  B bgf_;
  SU3 ptU_[allocORD];

 public:

  explicit BGptSU3(const B& bgf) : bgf_(bgf) { }
   
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
    for(int i = 0; i < PTORD; i++)  ptU_[i] *= z;
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
    for(int i = 0; i < PTORD; i++)  ptU_[i] /= z;
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
    for(int i = 0; i < PTORD; i++)  ptU_[i] += z[i];
    bgf_ += z.bgf();
    return *this;
  };
  BGptSU3& operator-=(const BGptSU3 &z) {
    for(int i = 0; i < PTORD; i++)  ptU_[i] -= z[i];
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
    SU3 tmp[allocORD];
    // Handle products that involve pert. orders > 1
    for (int i = 0; i < PTORD - 1; ++i)
      for (int j = 0; j < PTORD  - 1 - i; ++j)
	tmp[i + j + 1] += ptU_[j ] *A[i];
    // Include first order contributions
    for (int i = 0; i < PTORD; ++i)
      tmp[i] += bgf_.ApplyFromLeft(A[i]) +
	A.bgf_.ApplyFromRight(ptU_[i]);
    // Zeroth order.
    bgf_ *= A.bgf_;
    // copy
    //ptU_ = tmp;
    std::copy(tmp, tmp + allocORD, ptU_);
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
    for (int i = 0; i < PTORD; ++i)
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


  /*
  inline friend BGptSU3 log(const BGptSU3& U){
    BGptSU3 result = U, Bnew, res = U;
    double segno = -1, aux;
    result.flag = 0;

    for(int i = 2; i <= PTORD; i++){
     Bnew.zero();
     for(int iU = 1; iU < PTORD; iU++){
       for(int iB = i-1; iB <= PTORD - iU; iB++){
 	Bnew.ptU[iU+iB-1] += B.ptU[iB-1]*U.ptU[iU-1];
       }
     }
     for(int eq = 0; eq < PTORD; eq++){
       B.ptU[eq] = Bnew.ptU[eq];
     }
     aux = segno/(double)i;
     res += aux*B;
     segno = -segno;
    }
    res.flag = 0;
    
    return res;
  }
  
  inline friend BGptSU3 exp(const BGptSU3& A){
    BGptSU3 B = A, Bnew, res = A;
    double den;

    for(int i = 2; i <= PTORD; i++){
      
      Bnew.zero();
      for(int iA = 1; iA < PTORD; iA++){
        for(int iB = i-1; iB <= PTORD - iA; iB++){
	  Bnew.ptU[iA+iB-1] += B.ptU[iB-1]*A.ptU[iA-1];
        }
      }
      
      den = 1./(double)i;
      
      B = den*Bnew;
      res += B;
    }
    res.flag = 1;
    return res;
 }
  */

  //BGptSU3& reH();

  //void prout();

};


// --- end of ptSU3 declarations ---

#endif
