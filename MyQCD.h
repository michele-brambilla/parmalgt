#ifndef _MY_QCD_H
#define _MY_QCD_H

#include"MyMath.h"
#include"choices.h"

class Gluon;
class SpinColor;

class Gluon {

 public:
  SU3 U[dim];  

  Gluon(){};

  Gluon(const Gluon& A) {
    for (int i=0; i < dim; i++)
      U[i] = A.U[i];
  }

  int write(FILE *filept){
    if(fwrite(&U, SZ_DB*NC*NC*2, dim, filept)) return 0;
    return 1;
  }

  int read(FILE *filept){
    if(fread(&U, SZ_DB*9*2, dim, filept)) return 0;
    return 1;
  }

  Gluon& operator=(const Gluon& A) {
    for (int i=0; i < dim; i++){
      U[i] = A.U[i];
    }
    return *this;
  };

  Gluon operator-() {
    Gluon res;
    for (int i=0; i < dim; i++){
      res.U[i] = -U[i];
    }
    return res;
  };
  
  SU3 operator*(const Gluon&) const;

  SpinColor operator*(const SpinColor&) const;
  
  Gluon operator*(const Cplx&) const;

  void dag() {  for(int i = 0; i < dim; i++) U[i].dag(); }

  void prout(){
    for(int i = 0; i < dim; i++)
      U[i].prout();
  };

};

Gluon dag(const Gluon &U);


/************* end class gluon *****************/
/***********************************************/
/************* class spincolor *****************/


class SpinColor {

 public:
  CVector psi[dim];

  CVector& operator[] (const int& i) { return psi[i]; }
  const CVector& operator[] (const int& i) const { return psi[i]; }

  SpinColor(){};

  SpinColor(const SpinColor& P) {
    for (int i=0; i < dim; i++)
      psi[i] = P.psi[i];
  }

  int write(FILE *filept){
    if(fwrite(&psi, SZ_DB*3*2, dim, filept)) return 0;
    return 1;
  }

  int read(FILE *filept){
    if(fread(&psi, SZ_DB*3*2, dim, filept)) return 0;
    return 1;
  }

  void zeros() {
    memset(psi,0,dim*NC*sizeof(CVector));
/* #if (dim == 4 && NC == 3) */
/*     (psi[0].whr[0].re = psi[0].whr[0].im =  */
/*      psi[0].whr[1].re = psi[0].whr[1].im = */
/*      psi[0].whr[2].re = psi[0].whr[2].im = */
/*      psi[1].whr[0].re = psi[1].whr[0].im =  */
/*      psi[1].whr[1].re = psi[1].whr[1].im = */
/*      psi[1].whr[2].re = psi[1].whr[2].im = */
/*      psi[2].whr[0].re = psi[2].whr[0].im =  */
/*      psi[2].whr[1].re = psi[2].whr[1].im = */
/*      psi[2].whr[2].re = psi[2].whr[2].im = */
/*      psi[3].whr[0].re = psi[3].whr[0].im =  */
/*      psi[3].whr[1].re = psi[3].whr[1].im = */
/*      psi[3].whr[2].re = psi[3].whr[2].im = 0 ); */
/* #else */
/*       for (int mu = 0; mu < dim; mu++) */
/* 	for (int a = 0; a < NC; a++){ */
/* 	  psi[mu].whr[a].re = 0; */
/* 	  psi[mu].whr[a].im = 0; */
/* 	} */
/* #endif */
  }

  SpinColor& operator=(const SpinColor& P) {
#if dim == 4
    psi[0] = P.psi[0];
    psi[1] = P.psi[1];
    psi[2] = P.psi[2];
    psi[3] = P.psi[3];
#else
    for (int i=0; i < dim; i++)
      psi[i] = P.psi[i];
#endif
    return *this;
  }

  SpinColor operator-() const{
    SpinColor res;
#if dim == 4
    res.psi[0] = -psi[0];
    res.psi[1] = -psi[1];
    res.psi[2] = -psi[2];
    res.psi[3] = -psi[3];
#else
    for (int i=0; i < dim; i++)
      res.psi[i] = -psi[i];
#endif
    return res;
  }


  SpinColor operator+(const SpinColor& S) const{
    SpinColor sc;
#if dim == 4
    sc.psi[0] = psi[0] + S.psi[0];
    sc.psi[1] = psi[1] + S.psi[1];
    sc.psi[2] = psi[2] + S.psi[2];
    sc.psi[3] = psi[3] + S.psi[3];
#else
    for (int mu = 0; mu < dim; mu++) {
      sc.psi[mu] = psi[mu] + S.psi[mu];
    }
#endif
    return sc;
  }
  
  SpinColor operator-(const SpinColor& S) const{
    SpinColor sc; 
    for (int mu = 0; mu < dim; mu++) {
      sc.psi[mu] = psi[mu]- S.psi[mu];
    }
    return sc;
  }

  SpinColor& operator+=(const SpinColor& S) {
    for (int mu = 0; mu < dim; mu++) {
      psi[mu] += S.psi[mu];
    }
    return *this;
  }

  SpinColor& operator-=(const SpinColor& S) {
    for (int mu = 0; mu < dim; mu++) {
      psi[mu] -= S.psi[mu];
    }
    return *this;
  }
  
  Cplx operator*(const SpinColor&) const; // cosa restituisce in realta' questo?

  SpinColor operator*(const Cplx& z) const{
    SpinColor sc;
    for (int mu = 0; mu < dim; mu++) {
      sc.psi[mu] = z*psi[mu];
    }
    return sc;
  }

  SpinColor operator/(const Cplx& z) const{
    SpinColor sc;
    for (int mu = 0; mu < dim; mu++) {
      sc.psi[mu] = psi[mu]/z;
    }
    return sc;
  }

  SpinColor& operator*=(const Cplx& z) {
    for (int mu = 0; mu < dim; mu++) {
      psi[mu] *= z;
    }
    return *this;
  }

  SpinColor& operator/=(const Cplx& z) {
    for (int mu = 0; mu < dim; mu++) {
      psi[mu] /= z;
    }
    return *this;
  }

  SpinColor operator*(const double& x) const {
    SpinColor sc;
    for (int mu = 0; mu < dim; mu++) {
      sc.psi[mu] = psi[mu]*x;
    }
    return sc;
  }

  SpinColor operator/(const double& x) const {
    SpinColor sc;
    for (int mu = 0; mu < dim; mu++) {
      sc.psi[mu] = psi[mu]/x;
    }
    return sc;
  }

  SpinColor& operator*=(const double& x) {
    for (int mu = 0; mu < dim; mu++) {
      psi[mu] *= x;
    }
    return *this;
  }

  SpinColor& operator/=(const double& x) {
    for (int mu = 0; mu < dim; mu++) {
      psi[mu] /= x;
    }
    return *this;
  }

  SpinColor operator*(const SU3& U) const {
    SpinColor res;    
    for (int i = 0; i < dim; i++)
      res.psi[i] = psi[i] * U;
    return res;
  }


  friend SpinColor operator*(const Cplx& z, const SpinColor& S) {
    SpinColor sc;
    for (int mu = 0; mu < dim; mu++) {
      sc.psi[mu] = z * S.psi[mu];
    }
    return sc;
  }
  
  friend SpinColor operator*(const double& x, const SpinColor& S) {
    SpinColor sc;
    for (int mu = 0; mu < dim; mu++) {
      sc.psi[mu] = x * S.psi[mu];
    }
    return sc;
  }

  friend SpinColor operator*(const SU3& U, const SpinColor& S) {
    SpinColor res;    
    for (int i = 0; i < dim; i++)
      res.psi[i] = U * S.psi[i];
    return res;
  }


  void dag();

  void prout(){
    for(int i = 0; i < dim; i++)
      psi[i].prout();
  };

  void Gamma(int mu);

  SpinColor gmright(int mu);
  SpinColor gmleft(int mu) const;

  void uno_p_gmu(SpinColor&, int);
  void uno_m_gmu(SpinColor&, int);
};


SpinColor dag(const SpinColor &P);
SpinColor Gamma(int, const SpinColor &P);

// costruire operatori compositi!
// (bi- e trilineari)


#endif
