#ifndef _MY_QCD_H_
#define _MY_QCD_H_

#include"../MyMath.h"
#include"../choices.h"

class ptSU2;
class ptBoson;



class ptSU2{
 private:

 public:

  Cplx flag;    // 0 = matrice nell'algebra, 1 = vuoto banale, z3 ...  
  SU2 ptU[allocORD];

  ptSU2(Cplx i = 0){
    flag = i;
  }
  
  ~ptSU2(){};
 
  SU2* handle();

  int write(FILE *filept) {
    if(fwrite(&flag, SZ_DB*2, 1, filept) 
       && fwrite(&ptU, SZ_DB*2*2*2, allocORD, filept) ) return 0;
    return 1;
  }

  int read(FILE *filept) {
    if(fread(&flag, SZ_DB*2, 1, filept) 
       && fread(&ptU, SZ_DB*2*2*2, allocORD, filept) ) return 0;
    return 1;
  }

  void zero(){
    memset(ptU, 0, SZ_DB*2*2*2*PTORD);
    flag  = 0;
  };

  void id(){
    memset(ptU, 0, SZ_DB*2*NC*NC*PTORD);
    flag  = 1;
  };

  void Tr(Cplx *tt);  

  ptSU2& Trless();

  
  ptSU2& operator=(const ptSU2&);

  ptSU2 operator+(const ptSU2&) const ;

  ptSU2 operator-(const ptSU2&) const;

  ptSU2 operator*(const ptSU2& A) const{
    ptSU2 B;

    for(int i = 0; i < (PTORD-1); i++){
      for(int j = 0; j < (PTORD-1-i); j++){
	B.ptU[i+j+1] += ptU[j]*A.ptU[i];
      }
    }
  
    if( (flag == 0) && (A.flag == 0) ) {
      B.flag = 0;
      return B;
    }
    else {
      for(int i = 0; i < PTORD; i++){
	B.ptU[i] += A.flag*ptU[i] + flag*A.ptU[i];
      }
      
      B.flag = flag * A.flag;
    
      return B;
    }

  };
  

  ptSU2& operator+=(const ptSU2&);

  ptSU2& operator-=(const ptSU2&);

  ptSU2& operator*=(const ptSU2& A);

  ptSU2 operator*(const Cplx& z) const{
    ptSU2 B;
    for(int i = 0; i < PTORD; i++){
      B.ptU[i] = ptU[i]*z;
    }
    B.flag = z*flag;
    return (B);
  };
  
  ptSU2 operator/(const Cplx &z) const{
    ptSU2 B;
    for(int i = 0; i < PTORD; i++){
      B.ptU[i] = ptU[i]/z;
    }
    B.flag = flag/z;
    return B;
  };
  
  ptSU2& operator*=(const Cplx &z){
    for(int i = 0; i < PTORD; i++){
      ptU[i] *= z;
    }
    flag *= z;
    return *this;
  };
  
  ptSU2& operator/=(const Cplx &z){
    for(int i = 0; i < PTORD; i++){
      ptU[i] /= z;
    }
    flag /= z;
    return *this;
  };

  inline friend ptSU2 operator*(const Cplx&, const ptSU2&);

  ptSU2 operator*(const double&) const ;
  ptSU2 operator/(const double&) const ;

  ptSU2& operator*=(const double&);
  ptSU2& operator/=(const double&);

  inline friend ptSU2 operator*(const double& x, const ptSU2& U){
    ptSU2 B;
    for(int i = 0; i < PTORD; i++){
      B.ptU[i] = x*U.ptU[i];
    }
    B.flag = x*U.flag;
    return B;
  }

  friend ptSU2 dag(const ptSU2&);

  inline friend ptSU2 log(const ptSU2& U){
    ptSU2 B = U, Bnew, res = U;
    double segno = -1, aux;
    B.flag = 0;

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
  
  inline friend ptSU2 exp(const ptSU2& A){
    ptSU2 B = A, Bnew, res = A;
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

  ptSU2& reH();

  void prout();

};

// --- end of ptSU2 declarations ---




class ptBoson {

 public:
  ptSU2 U[dim];  

  ptBoson(){
    for (int mu = 0; mu < dim; mu++)
      U[mu].flag = 1;
  };
  
  ptBoson(const ptBoson& A);

  int write(FILE *filept) {
    for (int mu = 0; mu < dim; mu++)
      if (U[mu].write(filept)) return 0;
    return 1;
  }
  
  int read(FILE *filept) {
    for (int mu = 0; mu < dim; mu++)
      if (U[mu].read(filept)) return 0;
    return 1;
  }

  ptBoson& operator=(const ptBoson& A);
  
  ptSU2 operator*(const ptBoson& A) const;

/*   ptSpinColor operator*(const ptSpinColor& P) const; */
  
  void prout();

  friend ptBoson dag(const ptBoson& A);

};

// --- end of ptBoson declarations ---



/* class ptWVector { */

/* public: */

/*   CVector ptCV[PTORD + 1]; */

/*   int write(FILE *filept) { */
/*     if(fwrite(&ptCV, SZ_DB*2*NC, PTORD+1, filept) ) return 0; */
/*     return 1; */
/*   } */

/*   int read(FILE *filept) { */
/*     if(fread(&ptCV, SZ_DB*2*NC, PTORD+1, filept) ) return 0; */
/*     return 1; */
/*   } */
  
/*   ptCVector& operator=(const ptCVector &ptcv) { */
/*     for (int i = 0; i <= PTORD; i++) */
/*       ptCV[i] = ptcv.ptCV[i]; */
/*     return *this; */
/*   } */
  
/*   ptCVector operator+(const ptCVector &cv0) const{ */
/*     ptCVector cv1; */
/*     for (int i = 0; i <= PTORD; i++) */
/*       cv1.ptCV[i] = ptCV[i] + cv0.ptCV[i]; */
/*     return cv1; */
/*   } */
  
/*   ptCVector operator-(const ptCVector &cv0) const{ */
/*     ptCVector cv1; */
/*     for (int i = 0; i <= PTORD; i++) */
/*       cv1.ptCV[i] = ptCV[i] - cv0.ptCV[i]; */
/*     return cv1; */
/*   } */
  
/*   ptCVector& operator+=(const ptCVector &cv0){ */
/*     for (int i = 0; i <= PTORD; i++) */
/*       ptCV[i] += cv0.ptCV[i]; */
/*     return *this; */
/*   } */
  
/*   ptCVector& operator-=(const ptCVector &cv0) { */
/*     for (int i = 0; i <= PTORD; i++) */
/*       ptCV[i] += cv0.ptCV[i]; */
/*     return *this; */
/*   } */
  
/*   ptCVector operator*(const Cplx& z) const{ */
/*     ptCVector cv;     */
/*     for (int i = 0; i <= PTORD; i++) */
/*       cv.ptCV[i] = ptCV[i] * z; */
/*     return cv; */
/*   } */

/*   ptCVector operator/(const Cplx& z) const{ */
/*     ptCVector cv; */
/*     for (int i = 0; i <= PTORD; i++) */
/*       cv.ptCV[i] = ptCV[i] / z; */
/*     return cv; */
/*   } */

/*   ptCVector& operator*=(const Cplx& z) { */
/*     for (int i = 0; i < PTORD; i++) */
/*       ptCV[i] *= z; */
/*     return *this; */
/*   } */
  
/*   ptCVector& operator/=(const Cplx& z) { */
/*     for (int i = 0; i <= PTORD; i++) */
/*       ptCV[i] = ptCV[i] / z; */
/*     return *this; */
/*   } */
  
/*   friend ptCVector operator*(const Cplx& z, const ptCVector& cv0) { */
/*     ptCVector cv1; */
/*     for (int i = 0; i <= PTORD; i++) */
/*       cv1.ptCV[i] = cv0.ptCV[i] * z; */
/*     return cv1; */
/*   } */
  


/*   ptCVector operator*(ptSU2& A) const{ */
/*     ptCVector cv; */
/*     for(int i = 0; i < (PTORD-1); i++){ */
/*       for(int j = 0; j < (PTORD-1-i); j++){ */
/* 	cv.ptCV[i+j+1] += ptCV[j]*A.ptU[i]; */
/*       } */
/*     } */
/*     if(A.flag != 0) { */
/*       for(int i = 0; i < PTORD; i++){ */
/* 	cv.ptCV[i] += A.flag*ptCV[i]; */
/*       } */
/*     } */
/*     return cv; */
/*   } */

/*   ptCVector operator-() const{ */
/*     ptCVector cv1; */
/*     for (int i = 0; i <= PTORD; i++) */
/*       cv1.ptCV[i] = -ptCV[i]; */
/*     return cv1; */
/*   } */
  

/*   friend ptCVector operator*(ptSU2 &A, ptCVector& cv0){ */
/*     ptCVector cv1; */
/*     for(int i = 0; i < (PTORD-1); i++){ */
/*       for(int j = 0; j < (PTORD-1-i); j++){ */
/* 	cv1.ptCV[i+j+1] += A.ptU[j]*cv0.ptCV[i]; */
/*       } */
/*     } */
/*     if(A.flag != 0) { */
/*       for(int i = 0; i < PTORD; i++){ */
/* 	cv1.ptCV[i] += A.flag*cv0.ptCV[i]; */
/*       } */
/*     } */
    
/*     return cv1; */
/*   }; */
  
/* }; */


/* // ---- end ptCVector declarations ------- */



/* class ptSpinColor { */

/* public: */
/*   ptCVector psi[dim]; */
  
/*   ptSpinColor(){}; */

/*   int write(FILE *filept) { */
/*     for (int mu = 0; mu < dim; mu++) */
/*       if (psi[mu].write(filept)) return 0; */
/*     return 1; */
/*   } */
  
/*   int read(FILE *filept) { */
/*     for (int mu = 0; mu < dim; mu++) */
/*       if (psi[mu].read(filept)) return 0; */
/*     return 1; */
/*   } */
  
/*   ptSpinColor& operator=(const ptSpinColor &ptsc) { */
/* #if dim == 4 */
/*       psi[0] = ptsc.psi[0]; */
/*       psi[1] = ptsc.psi[1]; */
/*       psi[2] = ptsc.psi[2]; */
/*       psi[3] = ptsc.psi[3]; */
/* #else */
/*     for (int i = 0; i < dim; i++) */
/*       psi[i] = ptsc.psi[i]; */
/* #endif */
/*     return *this; */
/*   } */

/*   ptSpinColor operator+(const ptSpinColor &ptsc0) const{ */
/*     ptSpinColor ptsc1; */
/* #if dim == 4 */
/*     ptsc1.psi[0] = psi[0] + ptsc0.psi[0]; */
/*     ptsc1.psi[1] = psi[1] + ptsc0.psi[1]; */
/*     ptsc1.psi[2] = psi[2] + ptsc0.psi[2]; */
/*     ptsc1.psi[3] = psi[3] + ptsc0.psi[3]; */
/* #else */
/*     for (int i = 0; i < dim; i++) */
/*       ptsc1.psi[i] = psi[i] + ptsc0.psi[i]; */
/* #endif */
/*     return ptsc1; */
/*   } */

/*   ptSpinColor operator-(const ptSpinColor &ptsc0) const{ */
/*     ptSpinColor ptsc1; */
/* #if dim == 4 */
/*     ptsc1.psi[0] = psi[0] - ptsc0.psi[0]; */
/*     ptsc1.psi[1] = psi[1] - ptsc0.psi[1]; */
/*     ptsc1.psi[2] = psi[2] - ptsc0.psi[2]; */
/*     ptsc1.psi[3] = psi[3] - ptsc0.psi[3]; */
/* #else */
/*     for (int i = 0; i < dim; i++) */
/*       ptsc1.psi[i] = psi[i] - ptsc0.psi[i]; */
/* #endif */
/*     return ptsc1; */
/*   } */

/*   ptSpinColor operator*(const ptSpinColor& ptsc) const; */

/*   ptSpinColor& operator+=(const ptSpinColor &ptsc) { */
/* #if dim == 4 */
/*     psi[0] += ptsc.psi[0]; */
/*     psi[1] += ptsc.psi[1]; */
/*     psi[2] += ptsc.psi[2]; */
/*     psi[3] += ptsc.psi[3]; */
/* #else */
/*     for (int i = 0; i < dim; i++) */
/*       psi[i] += ptsc.psi[i]; */
/* #endif */
/*     return *this; */
/*   } */
  
/*   ptSpinColor& operator-=(const ptSpinColor &ptsc) { */
/* #if dim == 4 */
/*     psi[0] -= ptsc.psi[0]; */
/*     psi[1] -= ptsc.psi[1]; */
/*     psi[2] -= ptsc.psi[2]; */
/*     psi[3] -= ptsc.psi[3]; */
/* #else */
/*     for (int i = 0; i < dim; i++) */
/*       psi[i] -= ptsc.psi[i]; */
/* #endif */
/*     return *this; */
/*   } */
  
  
  
/*   ptSpinColor operator*(const Cplx& z) const{ */
/*     ptSpinColor ptsc; */
/* #if dim == 4 */
/*     ptsc.psi[0] = psi[0] * z; */
/*     ptsc.psi[1] = psi[1] * z; */
/*     ptsc.psi[2] = psi[2] * z; */
/*     ptsc.psi[3] = psi[3] * z; */
/* #else */
/*     for (int i = 0; i < dim; i++) */
/*       ptsc.psi[i] = psi[i] * z; */
/* #endif */
/*     return ptsc; */
/*   } */
  
/*   ptSpinColor operator/(const Cplx& z) const{ */
/*     ptSpinColor ptsc; */
/* #if dim == 4 */
/*     Cplx w = 1.0/z; */
/*     ptsc.psi[0] = psi[0] * w; */
/*     ptsc.psi[1] = psi[1] * w; */
/*     ptsc.psi[2] = psi[2] * w; */
/*     ptsc.psi[3] = psi[3] * w; */
/* #else */
/*     for (int i = 0; i < dim; i++) */
/*       ptsc.psi[i] = psi[i] / z; */
/* #endif */
/*     return ptsc; */
/*   } */
  
/*   ptSpinColor& operator*=(const Cplx& z) { */
/* #if dim == 4 */
/*     psi[0] = psi[0] * z; */
/*     psi[1] = psi[1] * z; */
/*     psi[2] = psi[2] * z; */
/*     psi[3] = psi[3] * z; */
/* #else */
/*     for (int i = 0; i < dim; i++) */
/*       psi[i] *= z; */
/* #endif */
/*     return *this; */
/*   } */
  
/*   ptSpinColor& operator/=(const Cplx& z) { */
/* #if dim == 4 */
/*     Cplx w = 1.0/z; */
/*     psi[0] = psi[0] * w; */
/*     psi[1] = psi[1] * w; */
/*     psi[2] = psi[2] * w; */
/*     psi[3] = psi[3] * w; */
/* #else */
/*     for (int i = 0; i < dim; i++) */
/*       psi[i] /= z; */
/* #endif */
/*     return *this; */
/*   } */
  
/*   friend ptSpinColor operator*(const Cplx& z, const ptSpinColor &ptsc0) { */
/*     ptSpinColor ptsc1; */
/* #if dim == 4 */
/*     ptsc1.psi[0] = ptsc0.psi[0] * z; */
/*     ptsc1.psi[1] = ptsc0.psi[1] * z; */
/*     ptsc1.psi[2] = ptsc0.psi[2] * z; */
/*     ptsc1.psi[3] = ptsc0.psi[3] * z; */
/* #else */
/*     for (int i = 0; i < dim; i++) */
/*       ptsc1.psi[i] = z * ptsc0.psi[i]; */
/* #endif */
/*     return ptsc1; */
/*   } */

/*   ptSpinColor operator-() const{ */
/*     ptSpinColor ptsc1; */
/* #if dim == 4 */
/*     ptsc1.psi[0] = -psi[0]; */
/*     ptsc1.psi[1] = -psi[1]; */
/*     ptsc1.psi[2] = -psi[2]; */
/*     ptsc1.psi[3] = -psi[3]; */
/* #else     */
/*     for (int i = 0; i < dim; i++) */
/*       ptsc1.psi[i] = -psi[i]; */
/* #endif */
/*     return ptsc1; */
/*   } */

/*   friend ptSpinColor dag(const ptSpinColor &ptsc); */

/*   SpinColor xSU2(const ptSU2& U, int ptord) { */
/*     SpinColor S; */

/*     for (int mu = 0; mu < dim; mu++) { */
/*       for (int i = 0; i < ptord; i++) { */
/* 	S.psi[mu] +=  psi[mu].ptCV[i] * U.ptU[ptord - i - 1]; */
/*       } */
/*     } */
/*     return S; */
/*   } */

/*   ptSpinColor gmleft(int mu) ; */
/*   ptSpinColor gmleft(int mu, int nu) ; */
  
/* };  */

/* // ****** end class ptSpinColor ********* */



#endif
