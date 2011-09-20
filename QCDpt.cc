#include<stdio.h>
#include"QCDpt.h"

#define SWAP(a,b) { tmp = a; a = b; b = tmp; }

int PTORD = allocORD;

SU3* ptSU3::handle(){
  return (SU3*)&ptU;
}

ptSU3& ptSU3::operator=(const ptSU3& A){
  flag  = A.flag;
    for(int i = 0; i < PTORD; i++){
      for(int j = 0; j < 9; j++){
        ptU[i].whr[j] = A.ptU[i].whr[j];
      }
    }
  return *this;
}


ptSU3 ptSU3::operator+(const ptSU3& A) const{
  ptSU3 B;
  for(int i = 0; i < PTORD; i++){
    B.ptU[i] = ptU[i] + A.ptU[i];
  }
  B.flag = flag + A.flag;

  return B;
}


ptSU3 ptSU3::operator-(const ptSU3& A) const{
  ptSU3 B;
  for(int i = 0; i < PTORD; i++){
    B.ptU[i] = ptU[i] - A.ptU[i];
  }
  
  B.flag = flag - A.flag;

  return B;
}

ptSU3& ptSU3::operator+=(const ptSU3& A) {

  for(int i = 0; i < PTORD; i++){
    ptU[i] += A.ptU[i];
  }

  flag += A.flag;
  return *this;
}


ptSU3& ptSU3::operator-=(const ptSU3& A){

  for(int i = 0; i < PTORD; i++){
    ptU[i] -= A.ptU[i];
  }
  flag -= A.flag;
  return *this;
}


ptSU3& ptSU3::operator*=(const ptSU3 &A){
  ptSU3 B;

  for(int i = 0; i < (PTORD-1); i++){
    for(int j = 0; j < (PTORD-1-i); j++){
      B.ptU[i+j+1] += ptU[j]*A.ptU[i];
    }
  }

  if( (flag == 0) && (A.flag == 0) ) {
    B.flag = 0;
    for(int i = 0; i < PTORD; i++){
      B.ptU[i] += A.flag*ptU[i] + flag*A.ptU[i];
    }
  }
  else{
    B.flag = flag * A.flag;
  }

  *this = B;
  return *this;
}


void ptSU3::Tr(Cplx *tt){
  for(int i = 0; i < PTORD; i++){
    tt[i+1] = (ptU[i].whr[0] +
	       ptU[i].whr[4] +
	       ptU[i].whr[8]);
  }
  tt[0] = 3.*flag;
}

ptSU3& ptSU3::Trless(){
  Cplx z;
  for(int i = 0; i < PTORD; i++){
    z = D3*(ptU[i].whr[0] +
	    ptU[i].whr[4] +
	    ptU[i].whr[8]);
    ptU[i].whr[0] -= z;
    ptU[i].whr[4] -= z;
    ptU[i].whr[8] -= z;
  }
  flag = 0;
  return *this;
}
  

ptSU3 operator*(const Cplx& z, const ptSU3 &U) {
  ptSU3 B;
  for(int i = 0; i < PTORD; i++){
    B.ptU[i] = z*U.ptU[i];
  }
  B.flag = z*U.flag;
  return B;
}


ptSU3 ptSU3::operator*(const double& x) const{
  ptSU3 B;
  for(int i = 0; i < PTORD; i++){
    B.ptU[i] = ptU[i]*x;
  }
  B.flag = x*flag;
  return B;
}

ptSU3 ptSU3::operator/(const double& x) const{
  ptSU3 B;
  for(int i = 0; i < PTORD; i++){
    B.ptU[i] = ptU[i]/x;
  }
  B.flag = flag/x;
  return B;
}


ptSU3& ptSU3::operator*=(const double& x){
  for(int i = 0; i < PTORD; i++){
    ptU[i] *= x;
  }
  flag *= x;
  return *this;
}

ptSU3& ptSU3::operator/=(const double& x){
  for(int i = 0; i < PTORD; i++){
    ptU[i] /= x;
  }
  flag /= x;
  return *this;
}


ptSU3 dag(const ptSU3& A){
  ptSU3 B;
  for(int i = 0; i < PTORD; i++){
    B.ptU[i] = dag(A.ptU[i]);
  }
  
  B.flag.re =  A.flag.re;
  B.flag.im = -A.flag.im;
  return B;
}


ptSU3& ptSU3::reH(){
  Cplx tr;
  for(int i = 0; i < PTORD; i++){
    ptU[i] -= dag(ptU[i]);
    ptU[i] = .5*ptU[i];
    tr = ptU[i].Tr()*D3;
    ptU[i].whr[0] -= tr;
    ptU[i].whr[4] -= tr;
    ptU[i].whr[8] -= tr;
  }
  flag = 0;
  return *this;
}


void ptSU3::prout(){
  flag.prout();
  printf("\n");
  for(int i = 0; i < PTORD; i++){
    ptU[i].prout();
  }
}


// ---- end of ptSU3 --------------- //




ptGluon& ptGluon::operator=(const ptGluon& A) {
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


ptGluon::ptGluon(const ptGluon& A) {
  for (int i=0; i < dim; i++)
    U[i] = A.U[i];
}



ptSU3 ptGluon::operator*(const ptGluon &A) const{
  ptSU3 res;
#if dim == 4
  res = U[0]*A.U[0] + U[1]*A.U[1] + U[2]*A.U[2] + U[3]*A.U[3] ;
#else
  for(int i = 0; i < dim; i++){
    res += U[i]*A.U[i];
  }
#endif
  return res;
}


void ptGluon::prout(){
  for(int i = 0; i < dim; i++)
    U[i].prout();
}


ptGluon dag(const ptGluon &A){  
  ptGluon B;
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


// ------------- end of ptGluon  -----------------//


ptSpinColor ptSpinColor::gmleft(int mu) {
  ptSpinColor sc;
  Cplx i(0,1), i_(0,-1);
  
  switch (mu) {
  case 0:
     sc.psi[0] = i_* psi[3];
     sc.psi[1] = i_* psi[2];
     sc.psi[2] = i * psi[1];
     sc.psi[3] = i * psi[0];
    break;
  case 1:
    sc.psi[0] = -psi[3];
    sc.psi[1] =  psi[2];
    sc.psi[2] =  psi[1];
    sc.psi[3] = -psi[0];
    break;
  case 2:
    sc.psi[0] = i_* psi[2];
    sc.psi[1] = i * psi[3];
    sc.psi[2] = i * psi[0];
    sc.psi[3] = i_* psi[1];
    break;
  case 3:
    sc.psi[0] = psi[2];
    sc.psi[1] = psi[3];
    sc.psi[2] = psi[0];
    sc.psi[3] = psi[1];
    break;
  case 5:
    sc.psi[0] =  psi[0];
    sc.psi[1] =  psi[1];
    sc.psi[2] = -psi[2];
    sc.psi[3] = -psi[3];
    break;
  }
  
  return sc;
}




// ptSpinColor ptSpinColor::gmleft(int mu) {
//   ptSpinColor sc;
//   Cplx i(0,1), i_(0,-1);
  
//   switch (mu) {
//   case 0:
//      sc.psi[0] =  psi[0];
//      sc.psi[1] =  psi[1];
//      sc.psi[2] = -psi[2];
//      sc.psi[3] = -psi[3];
//     break;
//   case 1:
//     sc.psi[0] = i_ * psi[3];
//     sc.psi[1] = i_ * psi[2];
//     sc.psi[2] = i  * psi[1];
//     sc.psi[3] = i  * psi[0];
//     break;
//   case 2:
//     sc.psi[0] = -psi[3];
//     sc.psi[1] =  psi[2];
//     sc.psi[2] =  psi[1];
//     sc.psi[3] = -psi[0];
//     break;
//   case 3:
//     sc.psi[0] = i_ * psi[2];
//     sc.psi[1] = i  * psi[3];
//     sc.psi[2] = i  * psi[0];
//     sc.psi[3] = i_ * psi[1];
//     break;
//   case 5:
//     sc.psi[0] = -psi[2];
//     sc.psi[1] = -psi[3];
//     sc.psi[2] = -psi[0];
//     sc.psi[3] = -psi[1];
//     break;
//   }
  
//   return sc;
// }


// ptSpinColor ptSpinColor::gmleft(int mu, int nu) {
//   ptSpinColor sc;
//   Cplx i(0,1), i_(0,-1);
  
//   switch (6*mu+nu) {
    
//   case 0:
//     return sc;
//     break;
    
//   case 1:
//     return -((*this).gmleft(nu,mu));
//     break;
    
//   case 2:
//     return -((*this).gmleft(nu,mu));
//     break;
    
//   case 3:
//     return -((*this).gmleft(nu,mu));
//     break;
    
//   case 5:
//     sc.psi[0] = -psi[2];
//     sc.psi[1] = -psi[3];
//     sc.psi[2] =  psi[0];
//     sc.psi[3] =  psi[1];
//     break;

//   case 6:
//     sc.psi[0] = i * psi[3];
//     sc.psi[1] = i * psi[2];
//     sc.psi[2] = i * psi[1];
//     sc.psi[3] = i * psi[0];
//     break;

//   case 7:
//     return sc;
//     break;

//   case 8:
//     sc.psi[0] = i  * psi[0];
//     sc.psi[1] = i_ * psi[1];
//     sc.psi[2] = i  * psi[2];
//     sc.psi[3] = i_ * psi[3];
//     break;

//   case 9:
//     sc.psi[0] = -psi[1];
//     sc.psi[1] =  psi[0];
//     sc.psi[2] = -psi[3];
//     sc.psi[3] =  psi[2];
//     break;

//   case 11:
//     sc.psi[0] = i  * psi[1];
//     sc.psi[1] = i  * psi[0];
//     sc.psi[2] = i_ * psi[3];
//     sc.psi[3] = i_ * psi[2];
//     break;

//   case 12:
//     sc.psi[0] =  psi[3];
//     sc.psi[1] = -psi[2];
//     sc.psi[2] =  psi[1];
//     sc.psi[3] = -psi[0];
//     break;

//   case 13:
//     return (-((*this).gmleft(nu,mu)));
//     break;

//   case 15:
//     sc.psi[0] = i * psi[1];
//     sc.psi[1] = i * psi[0];
//     sc.psi[2] = i * psi[3];
//     sc.psi[3] = i * psi[2];
//     break;

//   case 17:
//     sc.psi[0] =  psi[1];
//     sc.psi[1] = -psi[0];
//     sc.psi[2] = -psi[3];
//     sc.psi[3] =  psi[2];
//     break;

//   case 18:
//     sc.psi[0] = i  * psi[2];
//     sc.psi[1] = i_ * psi[3];
//     sc.psi[2] = i  * psi[0];
//     sc.psi[3] = i_ * psi[1];
//     break;

//   case 19:
//     return (-((*this).gmleft(nu,mu)));
//     break;

//   case 20:
//     return (-((*this).gmleft(nu,mu)));
//     break;

//   case 21:
//     return sc;
//     break;

//   case 23:
//     sc.psi[0] = i  * psi[0];
//     sc.psi[1] = i_ * psi[1];
//     sc.psi[2] = i_ * psi[2];
//     sc.psi[3] = i  * psi[3];
//     break;

//   default:
//     printf("Combinazione di matrici gamma non esistente.\n");
//     exit(1);

//   }
  
//   return sc;
// }



void ptSpinColor::uno_p_gmu(SpinColor& out, int mu, int ord){

  switch (mu) {
  case 0:
#ifdef __INTEL_INTRINSIC__
    out.psi[0].whr[0].m = _mm_addsub_pd(_mm_shuffle_pd(psi[0].whr[0].m,psi[0].whr[0].m,1), psi[3].whr[0].m );
    out.psi[0].whr[1].m = _mm_addsub_pd(_mm_shuffle_pd(psi[0].whr[1].m,psi[0].whr[1].m,1), psi[3].whr[1].m );
    out.psi[0].whr[2].m = _mm_addsub_pd(_mm_shuffle_pd(psi[0].whr[2].m,psi[0].whr[2].m,1), psi[3].whr[2].m );
    out.psi[0].whr[0].m = _mm_shuffle_pd(out.psi[0].whr[0].m,out.psi[0].whr[0].m,1);
    out.psi[0].whr[1].m = _mm_shuffle_pd(out.psi[0].whr[1].m,out.psi[0].whr[1].m,1);
    out.psi[0].whr[2].m = _mm_shuffle_pd(out.psi[0].whr[2].m,out.psi[0].whr[2].m,1);

    out.psi[1].whr[0].m = _mm_addsub_pd(_mm_shuffle_pd(psi[1].whr[0].m,psi[1].whr[0].m,1), psi[2].whr[0].m );
    out.psi[1].whr[1].m = _mm_addsub_pd(_mm_shuffle_pd(psi[1].whr[1].m,psi[1].whr[1].m,1), psi[2].whr[1].m );
    out.psi[1].whr[2].m = _mm_addsub_pd(_mm_shuffle_pd(psi[1].whr[2].m,psi[1].whr[2].m,1), psi[2].whr[2].m );
    out.psi[1].whr[0].m = _mm_shuffle_pd(out.psi[1].whr[0].m,out.psi[1].whr[0].m,1);
    out.psi[1].whr[1].m = _mm_shuffle_pd(out.psi[1].whr[1].m,out.psi[1].whr[1].m,1);
    out.psi[1].whr[2].m = _mm_shuffle_pd(out.psi[1].whr[2].m,out.psi[1].whr[2].m,1);

    out.psi[2].whr[0].m = _mm_addsub_pd(psi[2].whr[0].m, _mm_shuffle_pd(psi[1].whr[0].m,psi[1].whr[0].m,1));
    out.psi[2].whr[1].m = _mm_addsub_pd(psi[2].whr[1].m, _mm_shuffle_pd(psi[1].whr[1].m,psi[1].whr[1].m,1));
    out.psi[2].whr[2].m = _mm_addsub_pd(psi[2].whr[2].m, _mm_shuffle_pd(psi[1].whr[2].m,psi[1].whr[2].m,1));

    out.psi[3].whr[0].m = _mm_addsub_pd(psi[3].whr[0].m, _mm_shuffle_pd(psi[0].whr[0].m,psi[0].whr[0].m,1));
    out.psi[3].whr[1].m = _mm_addsub_pd(psi[3].whr[1].m, _mm_shuffle_pd(psi[0].whr[1].m,psi[0].whr[1].m,1));
    out.psi[3].whr[2].m = _mm_addsub_pd(psi[3].whr[2].m, _mm_shuffle_pd(psi[0].whr[2].m,psi[0].whr[2].m,1));
#else
    out.psi[0].whr[0].re = psi[0].ptCV[ord].whr[0].re + psi[3].ptCV[ord].whr[0].im;    out.psi[0].whr[0].im = psi[0].ptCV[ord].whr[0].im - psi[3].ptCV[ord].whr[0].re;
    out.psi[0].whr[1].re = psi[0].ptCV[ord].whr[1].re + psi[3].ptCV[ord].whr[1].im;    out.psi[0].whr[1].im = psi[0].ptCV[ord].whr[1].im - psi[3].ptCV[ord].whr[1].re;
    out.psi[0].whr[2].re = psi[0].ptCV[ord].whr[2].re + psi[3].ptCV[ord].whr[2].im;    out.psi[0].whr[2].im = psi[0].ptCV[ord].whr[2].im - psi[3].ptCV[ord].whr[2].re;
    out.psi[1].whr[0].re = psi[1].ptCV[ord].whr[0].re + psi[2].ptCV[ord].whr[0].im;    out.psi[1].whr[0].im = psi[1].ptCV[ord].whr[0].im - psi[2].ptCV[ord].whr[0].re;
    out.psi[1].whr[1].re = psi[1].ptCV[ord].whr[1].re + psi[2].ptCV[ord].whr[1].im;    out.psi[1].whr[1].im = psi[1].ptCV[ord].whr[1].im - psi[2].ptCV[ord].whr[1].re;
    out.psi[1].whr[2].re = psi[1].ptCV[ord].whr[2].re + psi[2].ptCV[ord].whr[2].im;    out.psi[1].whr[2].im = psi[1].ptCV[ord].whr[2].im - psi[2].ptCV[ord].whr[2].re;
    out.psi[2].whr[0].re = psi[2].ptCV[ord].whr[0].re - psi[1].ptCV[ord].whr[0].im;    out.psi[2].whr[0].im = psi[2].ptCV[ord].whr[0].im + psi[1].ptCV[ord].whr[0].re;
    out.psi[2].whr[1].re = psi[2].ptCV[ord].whr[1].re - psi[1].ptCV[ord].whr[1].im;    out.psi[2].whr[1].im = psi[2].ptCV[ord].whr[1].im + psi[1].ptCV[ord].whr[1].re;
    out.psi[2].whr[2].re = psi[2].ptCV[ord].whr[2].re - psi[1].ptCV[ord].whr[2].im;    out.psi[2].whr[2].im = psi[2].ptCV[ord].whr[2].im + psi[1].ptCV[ord].whr[2].re;
    out.psi[3].whr[0].re = psi[3].ptCV[ord].whr[0].re - psi[0].ptCV[ord].whr[0].im;    out.psi[3].whr[0].im = psi[3].ptCV[ord].whr[0].im + psi[0].ptCV[ord].whr[0].re;
    out.psi[3].whr[1].re = psi[3].ptCV[ord].whr[1].re - psi[0].ptCV[ord].whr[1].im;    out.psi[3].whr[1].im = psi[3].ptCV[ord].whr[1].im + psi[0].ptCV[ord].whr[1].re;
    out.psi[3].whr[2].re = psi[3].ptCV[ord].whr[2].re - psi[0].ptCV[ord].whr[2].im;    out.psi[3].whr[2].im = psi[3].ptCV[ord].whr[2].im + psi[0].ptCV[ord].whr[2].re;
#endif
    break;

  case 1:
#ifdef __INTEL_INTRINSIC__
    out.psi[0].whr[0].m = _mm_sub_pd(psi[0].whr[0].m, psi[3].whr[0].m);
    out.psi[0].whr[1].m = _mm_sub_pd(psi[0].whr[1].m, psi[3].whr[1].m);
    out.psi[0].whr[2].m = _mm_sub_pd(psi[0].whr[2].m, psi[3].whr[2].m);

    out.psi[1].whr[0].m = _mm_add_pd(psi[1].whr[0].m, psi[2].whr[0].m);
    out.psi[1].whr[1].m = _mm_add_pd(psi[1].whr[1].m, psi[2].whr[1].m);
    out.psi[1].whr[2].m = _mm_add_pd(psi[1].whr[2].m, psi[2].whr[2].m);

    out.psi[2].whr[0].m = _mm_add_pd(psi[2].whr[0].m, psi[1].whr[0].m);
    out.psi[2].whr[1].m = _mm_add_pd(psi[2].whr[1].m, psi[1].whr[1].m);
    out.psi[2].whr[2].m = _mm_add_pd(psi[2].whr[2].m, psi[1].whr[2].m);

    out.psi[3].whr[0].m = _mm_sub_pd(psi[3].whr[0].m, psi[0].whr[0].m);
    out.psi[3].whr[1].m = _mm_sub_pd(psi[3].whr[1].m, psi[0].whr[1].m);
    out.psi[3].whr[2].m = _mm_sub_pd(psi[3].whr[2].m, psi[0].whr[2].m);
#else
    out.psi[0].whr[0].re = psi[0].ptCV[ord].whr[0].re - psi[3].ptCV[ord].whr[0].re;    out.psi[0].whr[0].im = psi[0].ptCV[ord].whr[0].im - psi[3].ptCV[ord].whr[0].im;
    out.psi[0].whr[1].re = psi[0].ptCV[ord].whr[1].re - psi[3].ptCV[ord].whr[1].re;    out.psi[0].whr[1].im = psi[0].ptCV[ord].whr[1].im - psi[3].ptCV[ord].whr[1].im;
    out.psi[0].whr[2].re = psi[0].ptCV[ord].whr[2].re - psi[3].ptCV[ord].whr[2].re;    out.psi[0].whr[2].im = psi[0].ptCV[ord].whr[2].im - psi[3].ptCV[ord].whr[2].im;
    out.psi[1].whr[0].re = psi[1].ptCV[ord].whr[0].re + psi[2].ptCV[ord].whr[0].re;    out.psi[1].whr[0].im = psi[1].ptCV[ord].whr[0].im + psi[2].ptCV[ord].whr[0].im;
    out.psi[1].whr[1].re = psi[1].ptCV[ord].whr[1].re + psi[2].ptCV[ord].whr[1].re;    out.psi[1].whr[1].im = psi[1].ptCV[ord].whr[1].im + psi[2].ptCV[ord].whr[1].im;
    out.psi[1].whr[2].re = psi[1].ptCV[ord].whr[2].re + psi[2].ptCV[ord].whr[2].re;    out.psi[1].whr[2].im = psi[1].ptCV[ord].whr[2].im + psi[2].ptCV[ord].whr[2].im;
    out.psi[2].whr[0].re = psi[2].ptCV[ord].whr[0].re + psi[1].ptCV[ord].whr[0].re;    out.psi[2].whr[0].im = psi[2].ptCV[ord].whr[0].im + psi[1].ptCV[ord].whr[0].im;
    out.psi[2].whr[1].re = psi[2].ptCV[ord].whr[1].re + psi[1].ptCV[ord].whr[1].re;    out.psi[2].whr[1].im = psi[2].ptCV[ord].whr[1].im + psi[1].ptCV[ord].whr[1].im;
    out.psi[2].whr[2].re = psi[2].ptCV[ord].whr[2].re + psi[1].ptCV[ord].whr[2].re;    out.psi[2].whr[2].im = psi[2].ptCV[ord].whr[2].im + psi[1].ptCV[ord].whr[2].im;
    out.psi[3].whr[0].re = psi[3].ptCV[ord].whr[0].re - psi[0].ptCV[ord].whr[0].re;    out.psi[3].whr[0].im = psi[3].ptCV[ord].whr[0].im - psi[0].ptCV[ord].whr[0].im;
    out.psi[3].whr[1].re = psi[3].ptCV[ord].whr[1].re - psi[0].ptCV[ord].whr[1].re;    out.psi[3].whr[1].im = psi[3].ptCV[ord].whr[1].im - psi[0].ptCV[ord].whr[1].im;
    out.psi[3].whr[2].re = psi[3].ptCV[ord].whr[2].re - psi[0].ptCV[ord].whr[2].re;    out.psi[3].whr[2].im = psi[3].ptCV[ord].whr[2].im - psi[0].ptCV[ord].whr[2].im;
#endif
    break;

  case 2:
#ifdef __INTEL_INTRINSIC__
    out.psi[0].whr[0].m = _mm_addsub_pd(_mm_shuffle_pd(psi[0].whr[0].m,psi[0].whr[0].m,1), psi[2].whr[0].m );
    out.psi[0].whr[1].m = _mm_addsub_pd(_mm_shuffle_pd(psi[0].whr[1].m,psi[0].whr[1].m,1), psi[2].whr[1].m );
    out.psi[0].whr[2].m = _mm_addsub_pd(_mm_shuffle_pd(psi[0].whr[2].m,psi[0].whr[2].m,1), psi[2].whr[2].m );
    out.psi[0].whr[0].m = _mm_shuffle_pd(out.psi[0].whr[0].m,out.psi[0].whr[0].m,1);
    out.psi[0].whr[1].m = _mm_shuffle_pd(out.psi[0].whr[1].m,out.psi[0].whr[1].m,1);
    out.psi[0].whr[2].m = _mm_shuffle_pd(out.psi[0].whr[2].m,out.psi[0].whr[2].m,1);

    out.psi[1].whr[0].m = _mm_addsub_pd(psi[1].whr[0].m, _mm_shuffle_pd(psi[3].whr[0].m,psi[3].whr[0].m,1));
    out.psi[1].whr[1].m = _mm_addsub_pd(psi[1].whr[1].m, _mm_shuffle_pd(psi[3].whr[1].m,psi[3].whr[1].m,1));
    out.psi[1].whr[2].m = _mm_addsub_pd(psi[1].whr[2].m, _mm_shuffle_pd(psi[3].whr[2].m,psi[3].whr[2].m,1));

    out.psi[2].whr[0].m = _mm_addsub_pd(psi[2].whr[0].m, _mm_shuffle_pd(psi[0].whr[0].m,psi[0].whr[0].m,1));
    out.psi[2].whr[1].m = _mm_addsub_pd(psi[2].whr[1].m, _mm_shuffle_pd(psi[0].whr[1].m,psi[0].whr[1].m,1));
    out.psi[2].whr[2].m = _mm_addsub_pd(psi[2].whr[2].m, _mm_shuffle_pd(psi[0].whr[2].m,psi[0].whr[2].m,1));

    out.psi[3].whr[0].m = _mm_addsub_pd(_mm_shuffle_pd(psi[3].whr[0].m,psi[3].whr[0].m,1), psi[1].whr[0].m );
    out.psi[3].whr[1].m = _mm_addsub_pd(_mm_shuffle_pd(psi[3].whr[1].m,psi[3].whr[1].m,1), psi[1].whr[1].m );
    out.psi[3].whr[2].m = _mm_addsub_pd(_mm_shuffle_pd(psi[3].whr[2].m,psi[3].whr[2].m,1), psi[1].whr[2].m );
    out.psi[3].whr[0].m = _mm_shuffle_pd(out.psi[3].whr[0].m,out.psi[3].whr[0].m,1);
    out.psi[3].whr[1].m = _mm_shuffle_pd(out.psi[3].whr[1].m,out.psi[3].whr[1].m,1);
    out.psi[3].whr[2].m = _mm_shuffle_pd(out.psi[3].whr[2].m,out.psi[3].whr[2].m,1);
#else
    out.psi[0].whr[0].re = psi[0].ptCV[ord].whr[0].re + psi[2].ptCV[ord].whr[0].im;    out.psi[0].whr[0].im = psi[0].ptCV[ord].whr[0].im - psi[2].ptCV[ord].whr[0].re;
    out.psi[0].whr[1].re = psi[0].ptCV[ord].whr[1].re + psi[2].ptCV[ord].whr[1].im;    out.psi[0].whr[1].im = psi[0].ptCV[ord].whr[1].im - psi[2].ptCV[ord].whr[1].re;
    out.psi[0].whr[2].re = psi[0].ptCV[ord].whr[2].re + psi[2].ptCV[ord].whr[2].im;    out.psi[0].whr[2].im = psi[0].ptCV[ord].whr[2].im - psi[2].ptCV[ord].whr[2].re;
    out.psi[1].whr[0].re = psi[1].ptCV[ord].whr[0].re - psi[3].ptCV[ord].whr[0].im;    out.psi[1].whr[0].im = psi[1].ptCV[ord].whr[0].im + psi[3].ptCV[ord].whr[0].re;
    out.psi[1].whr[1].re = psi[1].ptCV[ord].whr[1].re - psi[3].ptCV[ord].whr[1].im;    out.psi[1].whr[1].im = psi[1].ptCV[ord].whr[1].im + psi[3].ptCV[ord].whr[1].re;
    out.psi[1].whr[2].re = psi[1].ptCV[ord].whr[2].re - psi[3].ptCV[ord].whr[2].im;    out.psi[1].whr[2].im = psi[1].ptCV[ord].whr[2].im + psi[3].ptCV[ord].whr[2].re;
    out.psi[2].whr[0].re = psi[2].ptCV[ord].whr[0].re - psi[0].ptCV[ord].whr[0].im;    out.psi[2].whr[0].im = psi[2].ptCV[ord].whr[0].im + psi[0].ptCV[ord].whr[0].re;
    out.psi[2].whr[1].re = psi[2].ptCV[ord].whr[1].re - psi[0].ptCV[ord].whr[1].im;    out.psi[2].whr[1].im = psi[2].ptCV[ord].whr[1].im + psi[0].ptCV[ord].whr[1].re;
    out.psi[2].whr[2].re = psi[2].ptCV[ord].whr[2].re - psi[0].ptCV[ord].whr[2].im;    out.psi[2].whr[2].im = psi[2].ptCV[ord].whr[2].im + psi[0].ptCV[ord].whr[2].re;
    out.psi[3].whr[0].re = psi[3].ptCV[ord].whr[0].re + psi[1].ptCV[ord].whr[0].im;    out.psi[3].whr[0].im = psi[3].ptCV[ord].whr[0].im - psi[1].ptCV[ord].whr[0].re;
    out.psi[3].whr[1].re = psi[3].ptCV[ord].whr[1].re + psi[1].ptCV[ord].whr[1].im;    out.psi[3].whr[1].im = psi[3].ptCV[ord].whr[1].im - psi[1].ptCV[ord].whr[1].re;
    out.psi[3].whr[2].re = psi[3].ptCV[ord].whr[2].re + psi[1].ptCV[ord].whr[2].im;    out.psi[3].whr[2].im = psi[3].ptCV[ord].whr[2].im - psi[1].ptCV[ord].whr[2].re;
#endif
    break;

  case 3:
#ifdef __INTEL_INTRINSIC__
    out.psi[0].whr[0].m = _mm_add_pd(psi[0].whr[0].m, psi[2].whr[0].m);
    out.psi[0].whr[1].m = _mm_add_pd(psi[0].whr[1].m, psi[2].whr[1].m);
    out.psi[0].whr[2].m = _mm_add_pd(psi[0].whr[2].m, psi[2].whr[2].m);
    out.psi[1].whr[0].m = _mm_add_pd(psi[1].whr[0].m, psi[3].whr[0].m);
    out.psi[1].whr[1].m = _mm_add_pd(psi[1].whr[1].m, psi[3].whr[1].m);
    out.psi[1].whr[2].m = _mm_add_pd(psi[1].whr[2].m, psi[3].whr[2].m);
    out.psi[2].whr[0].m = _mm_add_pd(psi[2].whr[0].m, psi[0].whr[0].m);
    out.psi[2].whr[1].m = _mm_add_pd(psi[2].whr[1].m, psi[0].whr[1].m);
    out.psi[2].whr[2].m = _mm_add_pd(psi[2].whr[2].m, psi[0].whr[2].m);
    out.psi[3].whr[0].m = _mm_add_pd(psi[3].whr[0].m, psi[1].whr[0].m);
    out.psi[3].whr[1].m = _mm_add_pd(psi[3].whr[1].m, psi[1].whr[1].m);
    out.psi[3].whr[2].m = _mm_add_pd(psi[3].whr[2].m, psi[1].whr[2].m);
#else
    out.psi[0].whr[0].re = psi[0].ptCV[ord].whr[0].re + psi[2].ptCV[ord].whr[0].re;    out.psi[0].whr[0].im = psi[0].ptCV[ord].whr[0].im + psi[2].ptCV[ord].whr[0].im;
    out.psi[0].whr[1].re = psi[0].ptCV[ord].whr[1].re + psi[2].ptCV[ord].whr[1].re;    out.psi[0].whr[1].im = psi[0].ptCV[ord].whr[1].im + psi[2].ptCV[ord].whr[1].im;
    out.psi[0].whr[2].re = psi[0].ptCV[ord].whr[2].re + psi[2].ptCV[ord].whr[2].re;    out.psi[0].whr[2].im = psi[0].ptCV[ord].whr[2].im + psi[2].ptCV[ord].whr[2].im;
    out.psi[1].whr[0].re = psi[1].ptCV[ord].whr[0].re + psi[3].ptCV[ord].whr[0].re;    out.psi[1].whr[0].im = psi[1].ptCV[ord].whr[0].im + psi[3].ptCV[ord].whr[0].im;
    out.psi[1].whr[1].re = psi[1].ptCV[ord].whr[1].re + psi[3].ptCV[ord].whr[1].re;    out.psi[1].whr[1].im = psi[1].ptCV[ord].whr[1].im + psi[3].ptCV[ord].whr[1].im;
    out.psi[1].whr[2].re = psi[1].ptCV[ord].whr[2].re + psi[3].ptCV[ord].whr[2].re;    out.psi[1].whr[2].im = psi[1].ptCV[ord].whr[2].im + psi[3].ptCV[ord].whr[2].im;
    out.psi[2].whr[0].re = psi[2].ptCV[ord].whr[0].re + psi[0].ptCV[ord].whr[0].re;    out.psi[2].whr[0].im = psi[2].ptCV[ord].whr[0].im + psi[0].ptCV[ord].whr[0].im;
    out.psi[2].whr[1].re = psi[2].ptCV[ord].whr[1].re + psi[0].ptCV[ord].whr[1].re;    out.psi[2].whr[1].im = psi[2].ptCV[ord].whr[1].im + psi[0].ptCV[ord].whr[1].im;
    out.psi[2].whr[2].re = psi[2].ptCV[ord].whr[2].re + psi[0].ptCV[ord].whr[2].re;    out.psi[2].whr[2].im = psi[2].ptCV[ord].whr[2].im + psi[0].ptCV[ord].whr[2].im;
    out.psi[3].whr[0].re = psi[3].ptCV[ord].whr[0].re + psi[1].ptCV[ord].whr[0].re;    out.psi[3].whr[0].im = psi[3].ptCV[ord].whr[0].im + psi[1].ptCV[ord].whr[0].im;
    out.psi[3].whr[1].re = psi[3].ptCV[ord].whr[1].re + psi[1].ptCV[ord].whr[1].re;    out.psi[3].whr[1].im = psi[3].ptCV[ord].whr[1].im + psi[1].ptCV[ord].whr[1].im;
    out.psi[3].whr[2].re = psi[3].ptCV[ord].whr[2].re + psi[1].ptCV[ord].whr[2].re;    out.psi[3].whr[2].im = psi[3].ptCV[ord].whr[2].im + psi[1].ptCV[ord].whr[2].im;
#endif
    break;
  }

} // uno_p_gmu


void ptSpinColor::uno_m_gmu(SpinColor& out, int mu, int ord){

  switch (mu) {
  case 0:
#ifdef __INTEL_INTRINSIC__
    out.psi[0].whr[0].m = _mm_addsub_pd(psi[0].whr[0].m, _mm_shuffle_pd(psi[3].whr[0].m, psi[3].whr[0].m,1));
    out.psi[0].whr[1].m = _mm_addsub_pd(psi[0].whr[1].m, _mm_shuffle_pd(psi[3].whr[1].m, psi[3].whr[1].m,1));
    out.psi[0].whr[2].m = _mm_addsub_pd(psi[0].whr[2].m, _mm_shuffle_pd(psi[3].whr[2].m, psi[3].whr[2].m,1));
    out.psi[1].whr[0].m = _mm_addsub_pd(psi[1].whr[0].m, _mm_shuffle_pd(psi[2].whr[0].m, psi[2].whr[0].m,1));
    out.psi[1].whr[1].m = _mm_addsub_pd(psi[1].whr[1].m, _mm_shuffle_pd(psi[2].whr[1].m, psi[2].whr[1].m,1));
    out.psi[1].whr[2].m = _mm_addsub_pd(psi[1].whr[2].m, _mm_shuffle_pd(psi[2].whr[2].m, psi[2].whr[2].m,1));

    out.psi[2].whr[0].m = _mm_addsub_pd(_mm_shuffle_pd(psi[2].whr[0].m,psi[2].whr[0].m,1), psi[1].whr[0].m );
    out.psi[2].whr[1].m = _mm_addsub_pd(_mm_shuffle_pd(psi[2].whr[1].m,psi[2].whr[1].m,1), psi[1].whr[1].m );
    out.psi[2].whr[2].m = _mm_addsub_pd(_mm_shuffle_pd(psi[2].whr[2].m,psi[2].whr[2].m,1), psi[1].whr[2].m );
    out.psi[2].whr[0].m = _mm_shuffle_pd(out.psi[2].whr[0].m,out.psi[2].whr[0].m,1);
    out.psi[2].whr[1].m = _mm_shuffle_pd(out.psi[2].whr[1].m,out.psi[2].whr[1].m,1);
    out.psi[2].whr[2].m = _mm_shuffle_pd(out.psi[2].whr[2].m,out.psi[2].whr[2].m,1);

    out.psi[3].whr[0].m = _mm_addsub_pd(_mm_shuffle_pd(psi[3].whr[0].m,psi[3].whr[0].m,1), psi[0].whr[0].m );
    out.psi[3].whr[1].m = _mm_addsub_pd(_mm_shuffle_pd(psi[3].whr[1].m,psi[3].whr[1].m,1), psi[0].whr[1].m );
    out.psi[3].whr[2].m = _mm_addsub_pd(_mm_shuffle_pd(psi[3].whr[2].m,psi[3].whr[2].m,1), psi[0].whr[2].m );
    out.psi[3].whr[0].m = _mm_shuffle_pd(out.psi[3].whr[0].m,out.psi[3].whr[0].m,1);
    out.psi[3].whr[1].m = _mm_shuffle_pd(out.psi[3].whr[1].m,out.psi[3].whr[1].m,1);
    out.psi[3].whr[2].m = _mm_shuffle_pd(out.psi[3].whr[2].m,out.psi[3].whr[2].m,1);
#else
    out.psi[0].whr[0].re = psi[0].ptCV[ord].whr[0].re - psi[3].ptCV[ord].whr[0].im;      out.psi[0].whr[0].im = psi[0].ptCV[ord].whr[0].im + psi[3].ptCV[ord].whr[0].re;  
    out.psi[0].whr[1].re = psi[0].ptCV[ord].whr[1].re - psi[3].ptCV[ord].whr[1].im;      out.psi[0].whr[1].im = psi[0].ptCV[ord].whr[1].im + psi[3].ptCV[ord].whr[1].re;  
    out.psi[0].whr[2].re = psi[0].ptCV[ord].whr[2].re - psi[3].ptCV[ord].whr[2].im;      out.psi[0].whr[2].im = psi[0].ptCV[ord].whr[2].im + psi[3].ptCV[ord].whr[2].re;  
    out.psi[1].whr[0].re = psi[1].ptCV[ord].whr[0].re - psi[2].ptCV[ord].whr[0].im;      out.psi[1].whr[0].im = psi[1].ptCV[ord].whr[0].im + psi[2].ptCV[ord].whr[0].re;  
    out.psi[1].whr[1].re = psi[1].ptCV[ord].whr[1].re - psi[2].ptCV[ord].whr[1].im;      out.psi[1].whr[1].im = psi[1].ptCV[ord].whr[1].im + psi[2].ptCV[ord].whr[1].re;  
    out.psi[1].whr[2].re = psi[1].ptCV[ord].whr[2].re - psi[2].ptCV[ord].whr[2].im;      out.psi[1].whr[2].im = psi[1].ptCV[ord].whr[2].im + psi[2].ptCV[ord].whr[2].re;  
    out.psi[2].whr[0].re = psi[2].ptCV[ord].whr[0].re + psi[1].ptCV[ord].whr[0].im;      out.psi[2].whr[0].im = psi[2].ptCV[ord].whr[0].im - psi[1].ptCV[ord].whr[0].re;  
    out.psi[2].whr[1].re = psi[2].ptCV[ord].whr[1].re + psi[1].ptCV[ord].whr[1].im;      out.psi[2].whr[1].im = psi[2].ptCV[ord].whr[1].im - psi[1].ptCV[ord].whr[1].re;  
    out.psi[2].whr[2].re = psi[2].ptCV[ord].whr[2].re + psi[1].ptCV[ord].whr[2].im;      out.psi[2].whr[2].im = psi[2].ptCV[ord].whr[2].im - psi[1].ptCV[ord].whr[2].re;  
    out.psi[3].whr[0].re = psi[3].ptCV[ord].whr[0].re + psi[0].ptCV[ord].whr[0].im;      out.psi[3].whr[0].im = psi[3].ptCV[ord].whr[0].im - psi[0].ptCV[ord].whr[0].re;  
    out.psi[3].whr[1].re = psi[3].ptCV[ord].whr[1].re + psi[0].ptCV[ord].whr[1].im;      out.psi[3].whr[1].im = psi[3].ptCV[ord].whr[1].im - psi[0].ptCV[ord].whr[1].re;  
    out.psi[3].whr[2].re = psi[3].ptCV[ord].whr[2].re + psi[0].ptCV[ord].whr[2].im;      out.psi[3].whr[2].im = psi[3].ptCV[ord].whr[2].im - psi[0].ptCV[ord].whr[2].re;  
#endif
    break;

  case 1:
#ifdef __INTEL_INTRINSIC__
    out.psi[0].whr[0].m = _mm_add_pd(psi[0].whr[0].m, psi[3].whr[0].m);
    out.psi[0].whr[1].m = _mm_add_pd(psi[0].whr[1].m, psi[3].whr[1].m);
    out.psi[0].whr[2].m = _mm_add_pd(psi[0].whr[2].m, psi[3].whr[2].m);
    out.psi[1].whr[0].m = _mm_sub_pd(psi[1].whr[0].m, psi[2].whr[0].m);
    out.psi[1].whr[1].m = _mm_sub_pd(psi[1].whr[1].m, psi[2].whr[1].m);
    out.psi[1].whr[2].m = _mm_sub_pd(psi[1].whr[2].m, psi[2].whr[2].m);
    out.psi[2].whr[0].m = _mm_sub_pd(psi[2].whr[0].m, psi[1].whr[0].m);
    out.psi[2].whr[1].m = _mm_sub_pd(psi[2].whr[1].m, psi[1].whr[1].m);
    out.psi[2].whr[2].m = _mm_sub_pd(psi[2].whr[2].m, psi[1].whr[2].m);
    out.psi[3].whr[0].m = _mm_add_pd(psi[3].whr[0].m, psi[0].whr[0].m);
    out.psi[3].whr[1].m = _mm_add_pd(psi[3].whr[1].m, psi[0].whr[1].m);
    out.psi[3].whr[2].m = _mm_add_pd(psi[3].whr[2].m, psi[0].whr[2].m);
#else
     out.psi[0].whr[0].re = psi[0].ptCV[ord].whr[0].re + psi[3].ptCV[ord].whr[0].re;    out.psi[0].whr[0].im = psi[0].ptCV[ord].whr[0].im + psi[3].ptCV[ord].whr[0].im;
     out.psi[0].whr[1].re = psi[0].ptCV[ord].whr[1].re + psi[3].ptCV[ord].whr[1].re;    out.psi[0].whr[1].im = psi[0].ptCV[ord].whr[1].im + psi[3].ptCV[ord].whr[1].im;
     out.psi[0].whr[2].re = psi[0].ptCV[ord].whr[2].re + psi[3].ptCV[ord].whr[2].re;    out.psi[0].whr[2].im = psi[0].ptCV[ord].whr[2].im + psi[3].ptCV[ord].whr[2].im;
     out.psi[1].whr[0].re = psi[1].ptCV[ord].whr[0].re - psi[2].ptCV[ord].whr[0].re;    out.psi[1].whr[0].im = psi[1].ptCV[ord].whr[0].im - psi[2].ptCV[ord].whr[0].im;
     out.psi[1].whr[1].re = psi[1].ptCV[ord].whr[1].re - psi[2].ptCV[ord].whr[1].re;    out.psi[1].whr[1].im = psi[1].ptCV[ord].whr[1].im - psi[2].ptCV[ord].whr[1].im;
     out.psi[1].whr[2].re = psi[1].ptCV[ord].whr[2].re - psi[2].ptCV[ord].whr[2].re;    out.psi[1].whr[2].im = psi[1].ptCV[ord].whr[2].im - psi[2].ptCV[ord].whr[2].im;
     out.psi[2].whr[0].re = psi[2].ptCV[ord].whr[0].re - psi[1].ptCV[ord].whr[0].re;    out.psi[2].whr[0].im = psi[2].ptCV[ord].whr[0].im - psi[1].ptCV[ord].whr[0].im;
     out.psi[2].whr[1].re = psi[2].ptCV[ord].whr[1].re - psi[1].ptCV[ord].whr[1].re;    out.psi[2].whr[1].im = psi[2].ptCV[ord].whr[1].im - psi[1].ptCV[ord].whr[1].im;
     out.psi[2].whr[2].re = psi[2].ptCV[ord].whr[2].re - psi[1].ptCV[ord].whr[2].re;    out.psi[2].whr[2].im = psi[2].ptCV[ord].whr[2].im - psi[1].ptCV[ord].whr[2].im;
     out.psi[3].whr[0].re = psi[3].ptCV[ord].whr[0].re + psi[0].ptCV[ord].whr[0].re;    out.psi[3].whr[0].im = psi[3].ptCV[ord].whr[0].im + psi[0].ptCV[ord].whr[0].im;
     out.psi[3].whr[1].re = psi[3].ptCV[ord].whr[1].re + psi[0].ptCV[ord].whr[1].re;    out.psi[3].whr[1].im = psi[3].ptCV[ord].whr[1].im + psi[0].ptCV[ord].whr[1].im;
     out.psi[3].whr[2].re = psi[3].ptCV[ord].whr[2].re + psi[0].ptCV[ord].whr[2].re;    out.psi[3].whr[2].im = psi[3].ptCV[ord].whr[2].im + psi[0].ptCV[ord].whr[2].im;
#endif
     break;

  case 2:
#ifdef __INTEL_INTRINSIC__
    out.psi[0].whr[0].m = _mm_addsub_pd(psi[0].whr[0].m, _mm_shuffle_pd(psi[2].whr[0].m, psi[2].whr[0].m,1));
    out.psi[0].whr[1].m = _mm_addsub_pd(psi[0].whr[1].m, _mm_shuffle_pd(psi[2].whr[1].m, psi[2].whr[1].m,1));
    out.psi[0].whr[2].m = _mm_addsub_pd(psi[0].whr[2].m, _mm_shuffle_pd(psi[2].whr[2].m, psi[2].whr[2].m,1));

    out.psi[1].whr[0].m = _mm_addsub_pd(_mm_shuffle_pd(psi[1].whr[0].m,psi[1].whr[0].m,1), psi[3].whr[0].m );
    out.psi[1].whr[1].m = _mm_addsub_pd(_mm_shuffle_pd(psi[1].whr[1].m,psi[1].whr[1].m,1), psi[3].whr[1].m );
    out.psi[1].whr[2].m = _mm_addsub_pd(_mm_shuffle_pd(psi[1].whr[2].m,psi[1].whr[2].m,1), psi[3].whr[2].m );
    out.psi[1].whr[0].m = _mm_shuffle_pd(out.psi[1].whr[0].m,out.psi[1].whr[0].m,1);
    out.psi[1].whr[1].m = _mm_shuffle_pd(out.psi[1].whr[1].m,out.psi[1].whr[1].m,1);
    out.psi[1].whr[2].m = _mm_shuffle_pd(out.psi[1].whr[2].m,out.psi[1].whr[2].m,1);

    out.psi[2].whr[0].m = _mm_addsub_pd(_mm_shuffle_pd(psi[2].whr[0].m,psi[2].whr[0].m,1), psi[0].whr[0].m );
    out.psi[2].whr[1].m = _mm_addsub_pd(_mm_shuffle_pd(psi[2].whr[1].m,psi[2].whr[1].m,1), psi[0].whr[1].m );
    out.psi[2].whr[2].m = _mm_addsub_pd(_mm_shuffle_pd(psi[2].whr[2].m,psi[2].whr[2].m,1), psi[0].whr[2].m );
    out.psi[2].whr[0].m = _mm_shuffle_pd(out.psi[2].whr[0].m,out.psi[2].whr[0].m,1);
    out.psi[2].whr[1].m = _mm_shuffle_pd(out.psi[2].whr[1].m,out.psi[2].whr[1].m,1);
    out.psi[2].whr[2].m = _mm_shuffle_pd(out.psi[2].whr[2].m,out.psi[2].whr[2].m,1);

    out.psi[3].whr[0].m = _mm_addsub_pd(psi[3].whr[0].m, _mm_shuffle_pd(psi[1].whr[0].m, psi[1].whr[0].m,1));
    out.psi[3].whr[1].m = _mm_addsub_pd(psi[3].whr[1].m, _mm_shuffle_pd(psi[1].whr[1].m, psi[1].whr[1].m,1));
    out.psi[3].whr[2].m = _mm_addsub_pd(psi[3].whr[2].m, _mm_shuffle_pd(psi[1].whr[2].m, psi[1].whr[2].m,1));
#else
    out.psi[0].whr[0].re = psi[0].ptCV[ord].whr[0].re - psi[2].ptCV[ord].whr[0].im;    out.psi[0].whr[0].im = psi[0].ptCV[ord].whr[0].im + psi[2].ptCV[ord].whr[0].re;
    out.psi[0].whr[1].re = psi[0].ptCV[ord].whr[1].re - psi[2].ptCV[ord].whr[1].im;    out.psi[0].whr[1].im = psi[0].ptCV[ord].whr[1].im + psi[2].ptCV[ord].whr[1].re;
    out.psi[0].whr[2].re = psi[0].ptCV[ord].whr[2].re - psi[2].ptCV[ord].whr[2].im;    out.psi[0].whr[2].im = psi[0].ptCV[ord].whr[2].im + psi[2].ptCV[ord].whr[2].re;
    out.psi[1].whr[0].re = psi[1].ptCV[ord].whr[0].re + psi[3].ptCV[ord].whr[0].im;    out.psi[1].whr[0].im = psi[1].ptCV[ord].whr[0].im - psi[3].ptCV[ord].whr[0].re;
    out.psi[1].whr[1].re = psi[1].ptCV[ord].whr[1].re + psi[3].ptCV[ord].whr[1].im;    out.psi[1].whr[1].im = psi[1].ptCV[ord].whr[1].im - psi[3].ptCV[ord].whr[1].re;
    out.psi[1].whr[2].re = psi[1].ptCV[ord].whr[2].re + psi[3].ptCV[ord].whr[2].im;    out.psi[1].whr[2].im = psi[1].ptCV[ord].whr[2].im - psi[3].ptCV[ord].whr[2].re;
    out.psi[2].whr[0].re = psi[2].ptCV[ord].whr[0].re + psi[0].ptCV[ord].whr[0].im;    out.psi[2].whr[0].im = psi[2].ptCV[ord].whr[0].im - psi[0].ptCV[ord].whr[0].re;
    out.psi[2].whr[1].re = psi[2].ptCV[ord].whr[1].re + psi[0].ptCV[ord].whr[1].im;    out.psi[2].whr[1].im = psi[2].ptCV[ord].whr[1].im - psi[0].ptCV[ord].whr[1].re;
    out.psi[2].whr[2].re = psi[2].ptCV[ord].whr[2].re + psi[0].ptCV[ord].whr[2].im;    out.psi[2].whr[2].im = psi[2].ptCV[ord].whr[2].im - psi[0].ptCV[ord].whr[2].re;
    out.psi[3].whr[0].re = psi[3].ptCV[ord].whr[0].re - psi[1].ptCV[ord].whr[0].im;    out.psi[3].whr[0].im = psi[3].ptCV[ord].whr[0].im + psi[1].ptCV[ord].whr[0].re;
    out.psi[3].whr[1].re = psi[3].ptCV[ord].whr[1].re - psi[1].ptCV[ord].whr[1].im;    out.psi[3].whr[1].im = psi[3].ptCV[ord].whr[1].im + psi[1].ptCV[ord].whr[1].re;
    out.psi[3].whr[2].re = psi[3].ptCV[ord].whr[2].re - psi[1].ptCV[ord].whr[2].im;    out.psi[3].whr[2].im = psi[3].ptCV[ord].whr[2].im + psi[1].ptCV[ord].whr[2].re;
#endif

    break;

  case 3:
#ifdef __INTEL_INTRINSIC__
    out.psi[0].whr[0].m = _mm_sub_pd(psi[0].whr[0].m, psi[2].whr[0].m);
    out.psi[0].whr[1].m = _mm_sub_pd(psi[0].whr[1].m, psi[2].whr[1].m);
    out.psi[0].whr[2].m = _mm_sub_pd(psi[0].whr[2].m, psi[2].whr[2].m);
    out.psi[1].whr[0].m = _mm_sub_pd(psi[1].whr[0].m, psi[3].whr[0].m);
    out.psi[1].whr[1].m = _mm_sub_pd(psi[1].whr[1].m, psi[3].whr[1].m);
    out.psi[1].whr[2].m = _mm_sub_pd(psi[1].whr[2].m, psi[3].whr[2].m);
    out.psi[2].whr[0].m = _mm_sub_pd(psi[2].whr[0].m, psi[0].whr[0].m);
    out.psi[2].whr[1].m = _mm_sub_pd(psi[2].whr[1].m, psi[0].whr[1].m);
    out.psi[2].whr[2].m = _mm_sub_pd(psi[2].whr[2].m, psi[0].whr[2].m);
    out.psi[3].whr[0].m = _mm_sub_pd(psi[3].whr[0].m, psi[1].whr[0].m);
    out.psi[3].whr[1].m = _mm_sub_pd(psi[3].whr[1].m, psi[1].whr[1].m);
    out.psi[3].whr[2].m = _mm_sub_pd(psi[3].whr[2].m, psi[1].whr[2].m);
#else
    out.psi[0].whr[0].re = psi[0].ptCV[ord].whr[0].re - psi[2].ptCV[ord].whr[0].re;    out.psi[0].whr[0].im = psi[0].ptCV[ord].whr[0].im - psi[2].ptCV[ord].whr[0].im;
    out.psi[0].whr[1].re = psi[0].ptCV[ord].whr[1].re - psi[2].ptCV[ord].whr[1].re;    out.psi[0].whr[1].im = psi[0].ptCV[ord].whr[1].im - psi[2].ptCV[ord].whr[1].im;
    out.psi[0].whr[2].re = psi[0].ptCV[ord].whr[2].re - psi[2].ptCV[ord].whr[2].re;    out.psi[0].whr[2].im = psi[0].ptCV[ord].whr[2].im - psi[2].ptCV[ord].whr[2].im;
    out.psi[1].whr[0].re = psi[1].ptCV[ord].whr[0].re - psi[3].ptCV[ord].whr[0].re;    out.psi[1].whr[0].im = psi[1].ptCV[ord].whr[0].im - psi[3].ptCV[ord].whr[0].im;
    out.psi[1].whr[1].re = psi[1].ptCV[ord].whr[1].re - psi[3].ptCV[ord].whr[1].re;    out.psi[1].whr[1].im = psi[1].ptCV[ord].whr[1].im - psi[3].ptCV[ord].whr[1].im;
    out.psi[1].whr[2].re = psi[1].ptCV[ord].whr[2].re - psi[3].ptCV[ord].whr[2].re;    out.psi[1].whr[2].im = psi[1].ptCV[ord].whr[2].im - psi[3].ptCV[ord].whr[2].im;
    out.psi[2].whr[0].re = psi[2].ptCV[ord].whr[0].re - psi[0].ptCV[ord].whr[0].re;    out.psi[2].whr[0].im = psi[2].ptCV[ord].whr[0].im - psi[0].ptCV[ord].whr[0].im;
    out.psi[2].whr[1].re = psi[2].ptCV[ord].whr[1].re - psi[0].ptCV[ord].whr[1].re;    out.psi[2].whr[1].im = psi[2].ptCV[ord].whr[1].im - psi[0].ptCV[ord].whr[1].im;
    out.psi[2].whr[2].re = psi[2].ptCV[ord].whr[2].re - psi[0].ptCV[ord].whr[2].re;    out.psi[2].whr[2].im = psi[2].ptCV[ord].whr[2].im - psi[0].ptCV[ord].whr[2].im;
    out.psi[3].whr[0].re = psi[3].ptCV[ord].whr[0].re - psi[1].ptCV[ord].whr[0].re;    out.psi[3].whr[0].im = psi[3].ptCV[ord].whr[0].im - psi[1].ptCV[ord].whr[0].im;
    out.psi[3].whr[1].re = psi[3].ptCV[ord].whr[1].re - psi[1].ptCV[ord].whr[1].re;    out.psi[3].whr[1].im = psi[3].ptCV[ord].whr[1].im - psi[1].ptCV[ord].whr[1].im;
    out.psi[3].whr[2].re = psi[3].ptCV[ord].whr[2].re - psi[1].ptCV[ord].whr[2].re;    out.psi[3].whr[2].im = psi[3].ptCV[ord].whr[2].im - psi[1].ptCV[ord].whr[2].im;
#endif
    break;
  }

}// uno_m_gmu
