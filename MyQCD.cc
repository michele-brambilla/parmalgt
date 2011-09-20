//#include<iostream>
#include"MyQCD.h"
#define SWAP(a,b) { CVector tmp = a; a = b; b = tmp; }

SU3 Gluon::operator*(const Gluon &A) const{
  SU3 res;
#if dim == 4
  res = U[0]*A.U[0] + U[1]*A.U[1] + U[2]*A.U[2] + U[3]*A.U[3] ;
#else
  for(int i = 0; i < dim; i++) res += U[i]*A.U[i];
#endif
  return res;
}


SpinColor Gluon::operator*(const SpinColor &P) const{
  SpinColor res;
#if dim == 4
  res.psi[0] = U[0]*P.psi[0] ;
  res.psi[1] = U[1]*P.psi[1] ;
  res.psi[2] = U[2]*P.psi[2] ;
  res.psi[3] = U[3]*P.psi[3] ;
#else
  for(int i = 0; i < dim; i++) res.psi[i] = U[i]*P.psi[i];
#endif
  return res;
}


Gluon Gluon::operator*(const Cplx& z) const{
  Gluon res;
#if dim == 4
  res.U[0] = U[0]*z;
  res.U[1] = U[1]*z;
  res.U[2] = U[2]*z;
  res.U[3] = U[3]*z;
#else
  for(int i = 0; i < dim; i++) res.psi[i] = U[i]*z;
#endif
  return res;
}


Gluon dag(const Gluon &A) {  
  Gluon res;
#if dim == 4
  res.U[0] = dag(A.U[0]); 
  res.U[1] = dag(A.U[1]); 
  res.U[2] = dag(A.U[2]); 
  res.U[3] = dag(A.U[3]); 
#else
  for(int i = 0; i < dim; i++) res.U[i] = dag(A.U[i]);
#endif
  return res;
}




void SpinColor::dag(){
  for(int i = 0; i < dim; i++) psi[i].dag();
}


SpinColor dag(const SpinColor &P){
  SpinColor res;
#if dim == 4
  res.psi[0] = dag(P.psi[0]);
  res.psi[1] = dag(P.psi[1]);
  res.psi[2] = dag(P.psi[2]);
  res.psi[3] = dag(P.psi[3]);
#else
  for(int i = 0; i < dim; i++) res.psi[i] = dag(P.psi[i]);
#endif
  return res;
}


Cplx SpinColor::operator*(const SpinColor &P) const{
  Cplx res;
#if dim == 4
  res = psi[0]*P.psi[0] + psi[1]*P.psi[1] + psi[2]*P.psi[2] + psi[3]*P.psi[3];
#else
  for(int i = 0; i < dim; i++)   res += psi[i]*P.psi[i];
#endif
  return res;
}


// ** Moltiplicazione a sinistra per gamma mu ** //

// Rappresentazione CHIRALE //

SpinColor SpinColor::gmleft(int mu) const{
  SpinColor sc;
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


void SpinColor::uno_p_gmu(SpinColor& out, int mu){

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
    out.psi[0].whr[0].re = psi[0].whr[0].re + psi[3].whr[0].im;    out.psi[0].whr[0].im = psi[0].whr[0].im - psi[3].whr[0].re;
    out.psi[0].whr[1].re = psi[0].whr[1].re + psi[3].whr[1].im;    out.psi[0].whr[1].im = psi[0].whr[1].im - psi[3].whr[1].re;
    out.psi[0].whr[2].re = psi[0].whr[2].re + psi[3].whr[2].im;    out.psi[0].whr[2].im = psi[0].whr[2].im - psi[3].whr[2].re;
    out.psi[1].whr[0].re = psi[1].whr[0].re + psi[2].whr[0].im;    out.psi[1].whr[0].im = psi[1].whr[0].im - psi[2].whr[0].re;
    out.psi[1].whr[1].re = psi[1].whr[1].re + psi[2].whr[1].im;    out.psi[1].whr[1].im = psi[1].whr[1].im - psi[2].whr[1].re;
    out.psi[1].whr[2].re = psi[1].whr[2].re + psi[2].whr[2].im;    out.psi[1].whr[2].im = psi[1].whr[2].im - psi[2].whr[2].re;
    out.psi[2].whr[0].re = psi[2].whr[0].re - psi[1].whr[0].im;    out.psi[2].whr[0].im = psi[2].whr[0].im + psi[1].whr[0].re;
    out.psi[2].whr[1].re = psi[2].whr[1].re - psi[1].whr[1].im;    out.psi[2].whr[1].im = psi[2].whr[1].im + psi[1].whr[1].re;
    out.psi[2].whr[2].re = psi[2].whr[2].re - psi[1].whr[2].im;    out.psi[2].whr[2].im = psi[2].whr[2].im + psi[1].whr[2].re;
    out.psi[3].whr[0].re = psi[3].whr[0].re - psi[0].whr[0].im;    out.psi[3].whr[0].im = psi[3].whr[0].im + psi[0].whr[0].re;
    out.psi[3].whr[1].re = psi[3].whr[1].re - psi[0].whr[1].im;    out.psi[3].whr[1].im = psi[3].whr[1].im + psi[0].whr[1].re;
    out.psi[3].whr[2].re = psi[3].whr[2].re - psi[0].whr[2].im;    out.psi[3].whr[2].im = psi[3].whr[2].im + psi[0].whr[2].re;
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
    out.psi[0].whr[0].re = psi[0].whr[0].re - psi[3].whr[0].re;    out.psi[0].whr[0].im = psi[0].whr[0].im - psi[3].whr[0].im;
    out.psi[0].whr[1].re = psi[0].whr[1].re - psi[3].whr[1].re;    out.psi[0].whr[1].im = psi[0].whr[1].im - psi[3].whr[1].im;
    out.psi[0].whr[2].re = psi[0].whr[2].re - psi[3].whr[2].re;    out.psi[0].whr[2].im = psi[0].whr[2].im - psi[3].whr[2].im;

    out.psi[1].whr[0].re = psi[1].whr[0].re + psi[2].whr[0].re;    out.psi[1].whr[0].im = psi[1].whr[0].im + psi[2].whr[0].im;
    out.psi[1].whr[1].re = psi[1].whr[1].re + psi[2].whr[1].re;    out.psi[1].whr[1].im = psi[1].whr[1].im + psi[2].whr[1].im;
    out.psi[1].whr[2].re = psi[1].whr[2].re + psi[2].whr[2].re;    out.psi[1].whr[2].im = psi[1].whr[2].im + psi[2].whr[2].im;

    out.psi[2].whr[0].re = psi[2].whr[0].re + psi[1].whr[0].re;    out.psi[2].whr[0].im = psi[2].whr[0].im + psi[1].whr[0].im;
    out.psi[2].whr[1].re = psi[2].whr[1].re + psi[1].whr[1].re;    out.psi[2].whr[1].im = psi[2].whr[1].im + psi[1].whr[1].im;
    out.psi[2].whr[2].re = psi[2].whr[2].re + psi[1].whr[2].re;    out.psi[2].whr[2].im = psi[2].whr[2].im + psi[1].whr[2].im;

    out.psi[3].whr[0].re = psi[3].whr[0].re - psi[0].whr[0].re;    out.psi[3].whr[0].im = psi[3].whr[0].im - psi[0].whr[0].im;
    out.psi[3].whr[1].re = psi[3].whr[1].re - psi[0].whr[1].re;    out.psi[3].whr[1].im = psi[3].whr[1].im - psi[0].whr[1].im;
    out.psi[3].whr[2].re = psi[3].whr[2].re - psi[0].whr[2].re;    out.psi[3].whr[2].im = psi[3].whr[2].im - psi[0].whr[2].im;
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
    out.psi[0].whr[0].re = psi[0].whr[0].re + psi[2].whr[0].im;    out.psi[0].whr[0].im = psi[0].whr[0].im - psi[2].whr[0].re;
    out.psi[0].whr[1].re = psi[0].whr[1].re + psi[2].whr[1].im;    out.psi[0].whr[1].im = psi[0].whr[1].im - psi[2].whr[1].re;
    out.psi[0].whr[2].re = psi[0].whr[2].re + psi[2].whr[2].im;    out.psi[0].whr[2].im = psi[0].whr[2].im - psi[2].whr[2].re;

    out.psi[1].whr[0].re = psi[1].whr[0].re - psi[3].whr[0].im;    out.psi[1].whr[0].im = psi[1].whr[0].im + psi[3].whr[0].re;
    out.psi[1].whr[1].re = psi[1].whr[1].re - psi[3].whr[1].im;    out.psi[1].whr[1].im = psi[1].whr[1].im + psi[3].whr[1].re;
    out.psi[1].whr[2].re = psi[1].whr[2].re - psi[3].whr[2].im;    out.psi[1].whr[2].im = psi[1].whr[2].im + psi[3].whr[2].re;

    out.psi[2].whr[0].re = psi[2].whr[0].re - psi[0].whr[0].im;    out.psi[2].whr[0].im = psi[2].whr[0].im + psi[0].whr[0].re;
    out.psi[2].whr[1].re = psi[2].whr[1].re - psi[0].whr[1].im;    out.psi[2].whr[1].im = psi[2].whr[1].im + psi[0].whr[1].re;
    out.psi[2].whr[2].re = psi[2].whr[2].re - psi[0].whr[2].im;    out.psi[2].whr[2].im = psi[2].whr[2].im + psi[0].whr[2].re;

    out.psi[3].whr[0].re = psi[3].whr[0].re + psi[1].whr[0].im;    out.psi[3].whr[0].im = psi[3].whr[0].im - psi[1].whr[0].re;
    out.psi[3].whr[1].re = psi[3].whr[1].re + psi[1].whr[1].im;    out.psi[3].whr[1].im = psi[3].whr[1].im - psi[1].whr[1].re;
    out.psi[3].whr[2].re = psi[3].whr[2].re + psi[1].whr[2].im;    out.psi[3].whr[2].im = psi[3].whr[2].im - psi[1].whr[2].re;
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
    out.psi[0].whr[0].re = psi[0].whr[0].re + psi[2].whr[0].re;    out.psi[0].whr[0].im = psi[0].whr[0].im + psi[2].whr[0].im;
    out.psi[0].whr[1].re = psi[0].whr[1].re + psi[2].whr[1].re;    out.psi[0].whr[1].im = psi[0].whr[1].im + psi[2].whr[1].im;
    out.psi[0].whr[2].re = psi[0].whr[2].re + psi[2].whr[2].re;    out.psi[0].whr[2].im = psi[0].whr[2].im + psi[2].whr[2].im;
    out.psi[1].whr[0].re = psi[1].whr[0].re + psi[3].whr[0].re;    out.psi[1].whr[0].im = psi[1].whr[0].im + psi[3].whr[0].im;
    out.psi[1].whr[1].re = psi[1].whr[1].re + psi[3].whr[1].re;    out.psi[1].whr[1].im = psi[1].whr[1].im + psi[3].whr[1].im;
    out.psi[1].whr[2].re = psi[1].whr[2].re + psi[3].whr[2].re;    out.psi[1].whr[2].im = psi[1].whr[2].im + psi[3].whr[2].im;
    out.psi[2].whr[0].re = psi[2].whr[0].re + psi[0].whr[0].re;    out.psi[2].whr[0].im = psi[2].whr[0].im + psi[0].whr[0].im;
    out.psi[2].whr[1].re = psi[2].whr[1].re + psi[0].whr[1].re;    out.psi[2].whr[1].im = psi[2].whr[1].im + psi[0].whr[1].im;
    out.psi[2].whr[2].re = psi[2].whr[2].re + psi[0].whr[2].re;    out.psi[2].whr[2].im = psi[2].whr[2].im + psi[0].whr[2].im;
    out.psi[3].whr[0].re = psi[3].whr[0].re + psi[1].whr[0].re;    out.psi[3].whr[0].im = psi[3].whr[0].im + psi[1].whr[0].im;
    out.psi[3].whr[1].re = psi[3].whr[1].re + psi[1].whr[1].re;    out.psi[3].whr[1].im = psi[3].whr[1].im + psi[1].whr[1].im;
    out.psi[3].whr[2].re = psi[3].whr[2].re + psi[1].whr[2].re;    out.psi[3].whr[2].im = psi[3].whr[2].im + psi[1].whr[2].im;
#endif
    break;
  }

} // uno_p_gmu


void SpinColor::uno_m_gmu(SpinColor& out, int mu){

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
    out.psi[0].whr[0].re = psi[0].whr[0].re - psi[3].whr[0].im;      out.psi[0].whr[0].im = psi[0].whr[0].im + psi[3].whr[0].re;  
    out.psi[0].whr[1].re = psi[0].whr[1].re - psi[3].whr[1].im;      out.psi[0].whr[1].im = psi[0].whr[1].im + psi[3].whr[1].re;  
    out.psi[0].whr[2].re = psi[0].whr[2].re - psi[3].whr[2].im;      out.psi[0].whr[2].im = psi[0].whr[2].im + psi[3].whr[2].re;  
    out.psi[1].whr[0].re = psi[1].whr[0].re - psi[2].whr[0].im;      out.psi[1].whr[0].im = psi[1].whr[0].im + psi[2].whr[0].re;  
    out.psi[1].whr[1].re = psi[1].whr[1].re - psi[2].whr[1].im;      out.psi[1].whr[1].im = psi[1].whr[1].im + psi[2].whr[1].re;  
    out.psi[1].whr[2].re = psi[1].whr[2].re - psi[2].whr[2].im;      out.psi[1].whr[2].im = psi[1].whr[2].im + psi[2].whr[2].re;  
											     
    out.psi[2].whr[0].re = psi[2].whr[0].re + psi[1].whr[0].im;      out.psi[2].whr[0].im = psi[2].whr[0].im - psi[1].whr[0].re;  
    out.psi[2].whr[1].re = psi[2].whr[1].re + psi[1].whr[1].im;      out.psi[2].whr[1].im = psi[2].whr[1].im - psi[1].whr[1].re;  
    out.psi[2].whr[2].re = psi[2].whr[2].re + psi[1].whr[2].im;      out.psi[2].whr[2].im = psi[2].whr[2].im - psi[1].whr[2].re;  
    out.psi[3].whr[0].re = psi[3].whr[0].re + psi[0].whr[0].im;      out.psi[3].whr[0].im = psi[3].whr[0].im - psi[0].whr[0].re;  
    out.psi[3].whr[1].re = psi[3].whr[1].re + psi[0].whr[1].im;      out.psi[3].whr[1].im = psi[3].whr[1].im - psi[0].whr[1].re;  
    out.psi[3].whr[2].re = psi[3].whr[2].re + psi[0].whr[2].im;      out.psi[3].whr[2].im = psi[3].whr[2].im - psi[0].whr[2].re;  
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
    out.psi[0].whr[0].re = psi[0].whr[0].re + psi[3].whr[0].re;    out.psi[0].whr[0].im = psi[0].whr[0].im + psi[3].whr[0].im;
    out.psi[0].whr[1].re = psi[0].whr[1].re + psi[3].whr[1].re;    out.psi[0].whr[1].im = psi[0].whr[1].im + psi[3].whr[1].im;
    out.psi[0].whr[2].re = psi[0].whr[2].re + psi[3].whr[2].re;    out.psi[0].whr[2].im = psi[0].whr[2].im + psi[3].whr[2].im;
    out.psi[1].whr[0].re = psi[1].whr[0].re - psi[2].whr[0].re;    out.psi[1].whr[0].im = psi[1].whr[0].im - psi[2].whr[0].im;
    out.psi[1].whr[1].re = psi[1].whr[1].re - psi[2].whr[1].re;    out.psi[1].whr[1].im = psi[1].whr[1].im - psi[2].whr[1].im;
    out.psi[1].whr[2].re = psi[1].whr[2].re - psi[2].whr[2].re;    out.psi[1].whr[2].im = psi[1].whr[2].im - psi[2].whr[2].im;
    out.psi[2].whr[0].re = psi[2].whr[0].re - psi[1].whr[0].re;    out.psi[2].whr[0].im = psi[2].whr[0].im - psi[1].whr[0].im;
    out.psi[2].whr[1].re = psi[2].whr[1].re - psi[1].whr[1].re;    out.psi[2].whr[1].im = psi[2].whr[1].im - psi[1].whr[1].im;
    out.psi[2].whr[2].re = psi[2].whr[2].re - psi[1].whr[2].re;    out.psi[2].whr[2].im = psi[2].whr[2].im - psi[1].whr[2].im;
    out.psi[3].whr[0].re = psi[3].whr[0].re + psi[0].whr[0].re;    out.psi[3].whr[0].im = psi[3].whr[0].im + psi[0].whr[0].im;
    out.psi[3].whr[1].re = psi[3].whr[1].re + psi[0].whr[1].re;    out.psi[3].whr[1].im = psi[3].whr[1].im + psi[0].whr[1].im;
    out.psi[3].whr[2].re = psi[3].whr[2].re + psi[0].whr[2].re;    out.psi[3].whr[2].im = psi[3].whr[2].im + psi[0].whr[2].im;
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
    out.psi[0].whr[0].re = psi[0].whr[0].re - psi[2].whr[0].im;    out.psi[0].whr[0].im = psi[0].whr[0].im + psi[2].whr[0].re;
    out.psi[0].whr[1].re = psi[0].whr[1].re - psi[2].whr[1].im;    out.psi[0].whr[1].im = psi[0].whr[1].im + psi[2].whr[1].re;
    out.psi[0].whr[2].re = psi[0].whr[2].re - psi[2].whr[2].im;    out.psi[0].whr[2].im = psi[0].whr[2].im + psi[2].whr[2].re;
    out.psi[1].whr[0].re = psi[1].whr[0].re + psi[3].whr[0].im;    out.psi[1].whr[0].im = psi[1].whr[0].im - psi[3].whr[0].re;
    out.psi[1].whr[1].re = psi[1].whr[1].re + psi[3].whr[1].im;    out.psi[1].whr[1].im = psi[1].whr[1].im - psi[3].whr[1].re;
    out.psi[1].whr[2].re = psi[1].whr[2].re + psi[3].whr[2].im;    out.psi[1].whr[2].im = psi[1].whr[2].im - psi[3].whr[2].re;
    out.psi[2].whr[0].re = psi[2].whr[0].re + psi[0].whr[0].im;    out.psi[2].whr[0].im = psi[2].whr[0].im - psi[0].whr[0].re;
    out.psi[2].whr[1].re = psi[2].whr[1].re + psi[0].whr[1].im;    out.psi[2].whr[1].im = psi[2].whr[1].im - psi[0].whr[1].re;
    out.psi[2].whr[2].re = psi[2].whr[2].re + psi[0].whr[2].im;    out.psi[2].whr[2].im = psi[2].whr[2].im - psi[0].whr[2].re;
    out.psi[3].whr[0].re = psi[3].whr[0].re - psi[1].whr[0].im;    out.psi[3].whr[0].im = psi[3].whr[0].im + psi[1].whr[0].re;
    out.psi[3].whr[1].re = psi[3].whr[1].re - psi[1].whr[1].im;    out.psi[3].whr[1].im = psi[3].whr[1].im + psi[1].whr[1].re;
    out.psi[3].whr[2].re = psi[3].whr[2].re - psi[1].whr[2].im;    out.psi[3].whr[2].im = psi[3].whr[2].im + psi[1].whr[2].re;
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
    out.psi[0].whr[0].re = psi[0].whr[0].re - psi[2].whr[0].re;    out.psi[0].whr[0].im = psi[0].whr[0].im - psi[2].whr[0].im;
    out.psi[0].whr[1].re = psi[0].whr[1].re - psi[2].whr[1].re;    out.psi[0].whr[1].im = psi[0].whr[1].im - psi[2].whr[1].im;
    out.psi[0].whr[2].re = psi[0].whr[2].re - psi[2].whr[2].re;    out.psi[0].whr[2].im = psi[0].whr[2].im - psi[2].whr[2].im;
    out.psi[1].whr[0].re = psi[1].whr[0].re - psi[3].whr[0].re;    out.psi[1].whr[0].im = psi[1].whr[0].im - psi[3].whr[0].im;
    out.psi[1].whr[1].re = psi[1].whr[1].re - psi[3].whr[1].re;    out.psi[1].whr[1].im = psi[1].whr[1].im - psi[3].whr[1].im;
    out.psi[1].whr[2].re = psi[1].whr[2].re - psi[3].whr[2].re;    out.psi[1].whr[2].im = psi[1].whr[2].im - psi[3].whr[2].im;
    out.psi[2].whr[0].re = psi[2].whr[0].re - psi[0].whr[0].re;    out.psi[2].whr[0].im = psi[2].whr[0].im - psi[0].whr[0].im;
    out.psi[2].whr[1].re = psi[2].whr[1].re - psi[0].whr[1].re;    out.psi[2].whr[1].im = psi[2].whr[1].im - psi[0].whr[1].im;
    out.psi[2].whr[2].re = psi[2].whr[2].re - psi[0].whr[2].re;    out.psi[2].whr[2].im = psi[2].whr[2].im - psi[0].whr[2].im;
    out.psi[3].whr[0].re = psi[3].whr[0].re - psi[1].whr[0].re;    out.psi[3].whr[0].im = psi[3].whr[0].im - psi[1].whr[0].im;
    out.psi[3].whr[1].re = psi[3].whr[1].re - psi[1].whr[1].re;    out.psi[3].whr[1].im = psi[3].whr[1].im - psi[1].whr[1].im;
    out.psi[3].whr[2].re = psi[3].whr[2].re - psi[1].whr[2].re;    out.psi[3].whr[2].im = psi[3].whr[2].im - psi[1].whr[2].im;
#endif
    break;
  }

} // uno_m_gmu
