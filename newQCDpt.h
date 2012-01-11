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

  template <class C> BGptSU3& operator*=(const C &z) {
    for(int i = 0; i < PTORD; i++)  ptU_[i] *= z;
    bgf_ *= z;
    return *this;
  };

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

  void randomize() {
    static MyRand r(1235431);
    for (int i = 0; i < PTORD; ++i)
      for (int j = 0; j < 3; ++j)
        for (int k = 0; k < 3; ++k)
          ptU_[i](j, k) = Cplx(r.Rand(), r.Rand());
  };

  template<class C> BGptSU3 operator*(const C& z) const {
    BGptSU3 result(*this);
    return result *= z;
  };
  
  template <class C> BGptSU3 operator/(const C& z) const{
    return *this * (1./z);
  };
  
  template <class C> BGptSU3& operator/=(const C &z){
    return *this *= 1./z;
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
