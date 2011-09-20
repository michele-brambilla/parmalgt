#include<stdio.h>
#include"SU2pt.h"

#define SWAP(a,b) { tmp = a; a = b; b = tmp; }

//int PTORD = allocORD;


SU2* ptSU2::handle(){
  return (SU2*)&ptU;
}

ptSU2& ptSU2::operator=(const ptSU2& A){
  flag  = A.flag;
    for(int i = 0; i < PTORD; i++){
      for(int j = 0; j < 4; j++){
        ptU[i].whr[j] = A.ptU[i].whr[j];
      }
    }
  return *this;
}


ptSU2 ptSU2::operator+(const ptSU2& A) const{
  ptSU2 B;
  for(int i = 0; i < PTORD; i++){
    B.ptU[i] = ptU[i] + A.ptU[i];
  }
  B.flag = flag + A.flag;

  return B;
}


ptSU2 ptSU2::operator-(const ptSU2& A) const{
  ptSU2 B;
  for(int i = 0; i < PTORD; i++){
    B.ptU[i] = ptU[i] - A.ptU[i];
  }
  
  B.flag = flag - A.flag;

  return B;
}

// ptSU2 ptSU2::operator*(ptSU2 A){
//   ptSU2 B;

//   for(int i = 0; i < (PTORD-1); i++){
//     for(int j = 0; j < (PTORD-1-i); j++){
//       B.ptU[i+j+1] += ptU[j]*A.ptU[i];
//     }
//   }

//   if(flag*A.flag != 0) {
//     for(int i = 0; i < PTORD; i++){
//       B.ptU[i] += A.flag*ptU[i] + flag*A.ptU[i];
//     }
//   }

//   B.flag = flag * A.flag;

//   return B;
// }



ptSU2& ptSU2::operator+=(const ptSU2& A) {

  for(int i = 0; i < PTORD; i++){
    ptU[i] += A.ptU[i];
  }

  flag += A.flag;
  return *this;
}


ptSU2& ptSU2::operator-=(const ptSU2& A){

  for(int i = 0; i < PTORD; i++){
    ptU[i] -= A.ptU[i];
  }
  flag -= A.flag;
  return *this;
}


ptSU2& ptSU2::operator*=(const ptSU2 &A){
  ptSU2 B;

  for(int i = 0; i < (PTORD-1); i++){
    for(int j = 0; j < (PTORD-1-i); j++){
      B.ptU[i+j+1] += ptU[j]*A.ptU[i];
    }
  }

  if(flag*A.flag != 0) {
    for(int i = 0; i < PTORD; i++){
      B.ptU[i] += A.flag*ptU[i] + flag*A.ptU[i];
    }
  }

  flag *= A.flag;
  *this = B;
  return *this;
}


void ptSU2::Tr(Cplx *tt){
  for(int i = 0; i < PTORD; i++){
    tt[i+1] = (ptU[i].whr[0] +
	       ptU[i].whr[3]);
  }
  tt[0] = 2.*flag;
}

ptSU2& ptSU2::Trless(){
  Cplx z;
  for(int i = 0; i < PTORD; i++){
    z = .5*(ptU[i].whr[0] +
	    ptU[i].whr[3]);
    ptU[i].whr[0] -= z;
    ptU[i].whr[3] -= z;
  }
  flag = 0;
  return *this;
}
  

ptSU2 operator*(const Cplx& z, const ptSU2 &U) {
  ptSU2 B;
  for(int i = 0; i < PTORD; i++){
    B.ptU[i] = z*U.ptU[i];
  }
  B.flag = z*U.flag;
  return B;
}


ptSU2 ptSU2::operator*(const double& x) const{
  ptSU2 B;
  for(int i = 0; i < PTORD; i++){
    B.ptU[i] = ptU[i]*x;
  }
  B.flag = x*flag;
  return B;
}

ptSU2 ptSU2::operator/(const double& x) const{
  ptSU2 B;
  for(int i = 0; i < PTORD; i++){
    B.ptU[i] = ptU[i]/x;
  }
  B.flag = flag/x;
  return B;
}


ptSU2& ptSU2::operator*=(const double& x){
  for(int i = 0; i < PTORD; i++){
    ptU[i] *= x;
  }
  flag *= x;
  return *this;
}

ptSU2& ptSU2::operator/=(const double& x){
  for(int i = 0; i < PTORD; i++){
    ptU[i] /= x;
  }
  flag /= x;
  return *this;
}


ptSU2 dag(const ptSU2& A){
  ptSU2 B;
  for(int i = 0; i < PTORD; i++){
    B.ptU[i] = dag(A.ptU[i]);
  }
  
  B.flag.re =  A.flag.re;
  B.flag.im = -A.flag.im;
  return B;
}


ptSU2& ptSU2::reH(){
  Cplx tr;
  for(int i = 0; i < PTORD; i++){
    ptU[i] -= dag(ptU[i]);
    ptU[i] = .5*ptU[i];
    tr = ptU[i].Tr()*.5;
    ptU[i].whr[0] -= tr;
    ptU[i].whr[3] -= tr;
  }
  flag = 0;
  return *this;
}


void ptSU2::prout(){
  flag.prout();
  printf("\n");
  for(int i = 0; i < PTORD; i++){
    ptU[i].prout();
  }
}


// ---- end of ptSU2 --------------- //




ptBoson& ptBoson::operator=(const ptBoson& A) {
#if dim == 4
  U[0] = A.U[0];
  U[1] = A.U[1];
  U[2] = A.U[2];
  U[3] = A.U[3];
#else
    for (int i=0; i < dim; i++)
      U[i] = A.U[i];
#endif
  return *this;
}


ptBoson::ptBoson(const ptBoson& A) {
  for (int i=0; i < dim; i++)
    U[i] = A.U[i];
}



ptSU2 ptBoson::operator*(const ptBoson &A) const{
  ptSU2 res;
#if dim == 4
  res = U[0]*A.U[0] + U[1]*A.U[1] + U[2]*A.U[2] + U[3]*A.U[3] ;
#else
  for(int i = 0; i < dim; i++){
    res += U[i]*A.U[i];
  }
#endif
  return res;
}


void ptBoson::prout(){
  for(int i = 0; i < dim; i++)
    U[i].prout();
}


ptBoson dag(const ptBoson &A){  
  ptBoson B;
#if dim == 4
    B.U[0] = dag(A.U[0]);
    B.U[1] = dag(A.U[1]);
    B.U[2] = dag(A.U[2]);
    B.U[3] = dag(A.U[3]);
#else
  for(int i = 0; i < dim; i++){
    B.U[i] = dag(A.U[i]);
  } 
#endif
  return B;
}


// ------------- end of ptBoson  -----------------//


// ptSpinColor ptSpinColor::gmleft(int mu) {
//   ptSpinColor sc;
//   Cplx i(0,1), i_(0,-1);
  
//   switch (mu) {
//   case 0:
//      sc.psi[0] = i_* psi[3];
//      sc.psi[1] = i_* psi[2];
//      sc.psi[2] = i * psi[1];
//      sc.psi[3] = i * psi[0];
//     break;
//   case 1:
//     sc.psi[0] = -psi[3];
//     sc.psi[1] =  psi[2];
//     sc.psi[2] =  psi[1];
//     sc.psi[3] = -psi[0];
//     break;
//   case 2:
//     sc.psi[0] = i_* psi[2];
//     sc.psi[1] = i * psi[3];
//     sc.psi[2] = i * psi[0];
//     sc.psi[3] = i_* psi[1];
//     break;
//   case 3:
//     sc.psi[0] = psi[2];
//     sc.psi[1] = psi[3];
//     sc.psi[2] = psi[0];
//     sc.psi[3] = psi[1];
//     break;
//   case 5:
//     sc.psi[0] =  psi[0];
//     sc.psi[1] =  psi[1];
//     sc.psi[2] = -psi[2];
//     sc.psi[3] = -psi[3];
//     break;
//   }
  
//   return sc;
// }




// // ptSpinColor ptSpinColor::gmleft(int mu) {
// //   ptSpinColor sc;
// //   Cplx i(0,1), i_(0,-1);
  
// //   switch (mu) {
// //   case 0:
// //      sc.psi[0] =  psi[0];
// //      sc.psi[1] =  psi[1];
// //      sc.psi[2] = -psi[2];
// //      sc.psi[3] = -psi[3];
// //     break;
// //   case 1:
// //     sc.psi[0] = i_ * psi[3];
// //     sc.psi[1] = i_ * psi[2];
// //     sc.psi[2] = i  * psi[1];
// //     sc.psi[3] = i  * psi[0];
// //     break;
// //   case 2:
// //     sc.psi[0] = -psi[3];
// //     sc.psi[1] =  psi[2];
// //     sc.psi[2] =  psi[1];
// //     sc.psi[3] = -psi[0];
// //     break;
// //   case 3:
// //     sc.psi[0] = i_ * psi[2];
// //     sc.psi[1] = i  * psi[3];
// //     sc.psi[2] = i  * psi[0];
// //     sc.psi[3] = i_ * psi[1];
// //     break;
// //   case 5:
// //     sc.psi[0] = -psi[2];
// //     sc.psi[1] = -psi[3];
// //     sc.psi[2] = -psi[0];
// //     sc.psi[3] = -psi[1];
// //     break;
// //   }
  
// //   return sc;
// // }


// // ptSpinColor ptSpinColor::gmleft(int mu, int nu) {
// //   ptSpinColor sc;
// //   Cplx i(0,1), i_(0,-1);
  
// //   switch (6*mu+nu) {
    
// //   case 0:
// //     return sc;
// //     break;
    
// //   case 1:
// //     return -((*this).gmleft(nu,mu));
// //     break;
    
// //   case 2:
// //     return -((*this).gmleft(nu,mu));
// //     break;
    
// //   case 3:
// //     return -((*this).gmleft(nu,mu));
// //     break;
    
// //   case 5:
// //     sc.psi[0] = -psi[2];
// //     sc.psi[1] = -psi[3];
// //     sc.psi[2] =  psi[0];
// //     sc.psi[3] =  psi[1];
// //     break;

// //   case 6:
// //     sc.psi[0] = i * psi[3];
// //     sc.psi[1] = i * psi[2];
// //     sc.psi[2] = i * psi[1];
// //     sc.psi[3] = i * psi[0];
// //     break;

// //   case 7:
// //     return sc;
// //     break;

// //   case 8:
// //     sc.psi[0] = i  * psi[0];
// //     sc.psi[1] = i_ * psi[1];
// //     sc.psi[2] = i  * psi[2];
// //     sc.psi[3] = i_ * psi[3];
// //     break;

// //   case 9:
// //     sc.psi[0] = -psi[1];
// //     sc.psi[1] =  psi[0];
// //     sc.psi[2] = -psi[3];
// //     sc.psi[3] =  psi[2];
// //     break;

// //   case 11:
// //     sc.psi[0] = i  * psi[1];
// //     sc.psi[1] = i  * psi[0];
// //     sc.psi[2] = i_ * psi[3];
// //     sc.psi[3] = i_ * psi[2];
// //     break;

// //   case 12:
// //     sc.psi[0] =  psi[3];
// //     sc.psi[1] = -psi[2];
// //     sc.psi[2] =  psi[1];
// //     sc.psi[3] = -psi[0];
// //     break;

// //   case 13:
// //     return (-((*this).gmleft(nu,mu)));
// //     break;

// //   case 15:
// //     sc.psi[0] = i * psi[1];
// //     sc.psi[1] = i * psi[0];
// //     sc.psi[2] = i * psi[3];
// //     sc.psi[3] = i * psi[2];
// //     break;

// //   case 17:
// //     sc.psi[0] =  psi[1];
// //     sc.psi[1] = -psi[0];
// //     sc.psi[2] = -psi[3];
// //     sc.psi[3] =  psi[2];
// //     break;

// //   case 18:
// //     sc.psi[0] = i  * psi[2];
// //     sc.psi[1] = i_ * psi[3];
// //     sc.psi[2] = i  * psi[0];
// //     sc.psi[3] = i_ * psi[1];
// //     break;

// //   case 19:
// //     return (-((*this).gmleft(nu,mu)));
// //     break;

// //   case 20:
// //     return (-((*this).gmleft(nu,mu)));
// //     break;

// //   case 21:
// //     return sc;
// //     break;

// //   case 23:
// //     sc.psi[0] = i  * psi[0];
// //     sc.psi[1] = i_ * psi[1];
// //     sc.psi[2] = i_ * psi[2];
// //     sc.psi[3] = i  * psi[3];
// //     break;

// //   default:
// //     printf("Combinazione di matrici gamma non esistente.\n");
// //     exit(1);

// //   }
  
// //   return sc;
// // }
