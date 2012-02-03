#include <newQCDpt.h>

template  <int ORD, int DIM>
void BGptSpinColor<ORD, DIM>::uno_p_gmu(SpinColor& out, int mu, int ord){

  switch (mu) {
  case 0:
#ifdef __INTEL_INTRINSIC__
    out[0][0].m = _mm_addsub_pd(_mm_shuffle_pd(psi_[0][0].m,psi_[0][0].m,1), psi_[3][0].m );
    out[0][1].m = _mm_addsub_pd(_mm_shuffle_pd(psi_[0][1].m,psi_[0][1].m,1), psi_[3][1].m );
    out[0][2].m = _mm_addsub_pd(_mm_shuffle_pd(psi_[0][2].m,psi_[0][2].m,1), psi_[3][2].m );
    out[0][0].m = _mm_shuffle_pd(out[0][0].m,out[0][0].m,1);
    out[0][1].m = _mm_shuffle_pd(out[0][1].m,out[0][1].m,1);
    out[0][2].m = _mm_shuffle_pd(out[0][2].m,out[0][2].m,1);

    out[1][0].m = _mm_addsub_pd(_mm_shuffle_pd(psi_[1][0].m,psi_[1][0].m,1), psi_[2][0].m );
    out[1][1].m = _mm_addsub_pd(_mm_shuffle_pd(psi_[1][1].m,psi_[1][1].m,1), psi_[2][1].m );
    out[1][2].m = _mm_addsub_pd(_mm_shuffle_pd(psi_[1][2].m,psi_[1][2].m,1), psi_[2][2].m );
    out[1][0].m = _mm_shuffle_pd(out[1][0].m,out[1][0].m,1);
    out[1][1].m = _mm_shuffle_pd(out[1][1].m,out[1][1].m,1);
    out[1][2].m = _mm_shuffle_pd(out[1][2].m,out[1][2].m,1);

    out[2][0].m = _mm_addsub_pd(psi_[2][0].m, _mm_shuffle_pd(psi_[1][0].m,psi_[1][0].m,1));
    out[2][1].m = _mm_addsub_pd(psi_[2][1].m, _mm_shuffle_pd(psi_[1][1].m,psi_[1][1].m,1));
    out[2][2].m = _mm_addsub_pd(psi_[2][2].m, _mm_shuffle_pd(psi_[1][2].m,psi_[1][2].m,1));

    out[3][0].m = _mm_addsub_pd(psi_[3][0].m, _mm_shuffle_pd(psi_[0][0].m,psi_[0][0].m,1));
    out[3][1].m = _mm_addsub_pd(psi_[3][1].m, _mm_shuffle_pd(psi_[0][1].m,psi_[0][1].m,1));
    out[3][2].m = _mm_addsub_pd(psi_[3][2].m, _mm_shuffle_pd(psi_[0][2].m,psi_[0][2].m,1));
#else
    out[0][0].re = psi_[0][ord][0].re + psi_[3][ord][0].im;    out[0][0].im = psi_[0][ord][0].im - psi_[3][ord][0].re;
    out[0][1].re = psi_[0][ord][1].re + psi_[3][ord][1].im;    out[0][1].im = psi_[0][ord][1].im - psi_[3][ord][1].re;
    out[0][2].re = psi_[0][ord][2].re + psi_[3][ord][2].im;    out[0][2].im = psi_[0][ord][2].im - psi_[3][ord][2].re;
    out[1][0].re = psi_[1][ord][0].re + psi_[2][ord][0].im;    out[1][0].im = psi_[1][ord][0].im - psi_[2][ord][0].re;
    out[1][1].re = psi_[1][ord][1].re + psi_[2][ord][1].im;    out[1][1].im = psi_[1][ord][1].im - psi_[2][ord][1].re;
    out[1][2].re = psi_[1][ord][2].re + psi_[2][ord][2].im;    out[1][2].im = psi_[1][ord][2].im - psi_[2][ord][2].re;
    out[2][0].re = psi_[2][ord][0].re - psi_[1][ord][0].im;    out[2][0].im = psi_[2][ord][0].im + psi_[1][ord][0].re;
    out[2][1].re = psi_[2][ord][1].re - psi_[1][ord][1].im;    out[2][1].im = psi_[2][ord][1].im + psi_[1][ord][1].re;
    out[2][2].re = psi_[2][ord][2].re - psi_[1][ord][2].im;    out[2][2].im = psi_[2][ord][2].im + psi_[1][ord][2].re;
    out[3][0].re = psi_[3][ord][0].re - psi_[0][ord][0].im;    out[3][0].im = psi_[3][ord][0].im + psi_[0][ord][0].re;
    out[3][1].re = psi_[3][ord][1].re - psi_[0][ord][1].im;    out[3][1].im = psi_[3][ord][1].im + psi_[0][ord][1].re;
    out[3][2].re = psi_[3][ord][2].re - psi_[0][ord][2].im;    out[3][2].im = psi_[3][ord][2].im + psi_[0][ord][2].re;
#endif
    break;

  case 1:
#ifdef __INTEL_INTRINSIC__
    out[0][0].m = _mm_sub_pd(psi_[0][0].m, psi_[3][0].m);
    out[0][1].m = _mm_sub_pd(psi_[0][1].m, psi_[3][1].m);
    out[0][2].m = _mm_sub_pd(psi_[0][2].m, psi_[3][2].m);

    out[1][0].m = _mm_add_pd(psi_[1][0].m, psi_[2][0].m);
    out[1][1].m = _mm_add_pd(psi_[1][1].m, psi_[2][1].m);
    out[1][2].m = _mm_add_pd(psi_[1][2].m, psi_[2][2].m);

    out[2][0].m = _mm_add_pd(psi_[2][0].m, psi_[1][0].m);
    out[2][1].m = _mm_add_pd(psi_[2][1].m, psi_[1][1].m);
    out[2][2].m = _mm_add_pd(psi_[2][2].m, psi_[1][2].m);

    out[3][0].m = _mm_sub_pd(psi_[3][0].m, psi_[0][0].m);
    out[3][1].m = _mm_sub_pd(psi_[3][1].m, psi_[0][1].m);
    out[3][2].m = _mm_sub_pd(psi_[3][2].m, psi_[0][2].m);
#else
    out[0][0].re = psi_[0][ord][0].re - psi_[3][ord][0].re;    out[0][0].im = psi_[0][ord][0].im - psi_[3][ord][0].im;
    out[0][1].re = psi_[0][ord][1].re - psi_[3][ord][1].re;    out[0][1].im = psi_[0][ord][1].im - psi_[3][ord][1].im;
    out[0][2].re = psi_[0][ord][2].re - psi_[3][ord][2].re;    out[0][2].im = psi_[0][ord][2].im - psi_[3][ord][2].im;
    out[1][0].re = psi_[1][ord][0].re + psi_[2][ord][0].re;    out[1][0].im = psi_[1][ord][0].im + psi_[2][ord][0].im;
    out[1][1].re = psi_[1][ord][1].re + psi_[2][ord][1].re;    out[1][1].im = psi_[1][ord][1].im + psi_[2][ord][1].im;
    out[1][2].re = psi_[1][ord][2].re + psi_[2][ord][2].re;    out[1][2].im = psi_[1][ord][2].im + psi_[2][ord][2].im;
    out[2][0].re = psi_[2][ord][0].re + psi_[1][ord][0].re;    out[2][0].im = psi_[2][ord][0].im + psi_[1][ord][0].im;
    out[2][1].re = psi_[2][ord][1].re + psi_[1][ord][1].re;    out[2][1].im = psi_[2][ord][1].im + psi_[1][ord][1].im;
    out[2][2].re = psi_[2][ord][2].re + psi_[1][ord][2].re;    out[2][2].im = psi_[2][ord][2].im + psi_[1][ord][2].im;
    out[3][0].re = psi_[3][ord][0].re - psi_[0][ord][0].re;    out[3][0].im = psi_[3][ord][0].im - psi_[0][ord][0].im;
    out[3][1].re = psi_[3][ord][1].re - psi_[0][ord][1].re;    out[3][1].im = psi_[3][ord][1].im - psi_[0][ord][1].im;
    out[3][2].re = psi_[3][ord][2].re - psi_[0][ord][2].re;    out[3][2].im = psi_[3][ord][2].im - psi_[0][ord][2].im;
#endif
    break;

  case 2:
#ifdef __INTEL_INTRINSIC__
    out[0][0].m = _mm_addsub_pd(_mm_shuffle_pd(psi_[0][0].m,psi_[0][0].m,1), psi_[2][0].m );
    out[0][1].m = _mm_addsub_pd(_mm_shuffle_pd(psi_[0][1].m,psi_[0][1].m,1), psi_[2][1].m );
    out[0][2].m = _mm_addsub_pd(_mm_shuffle_pd(psi_[0][2].m,psi_[0][2].m,1), psi_[2][2].m );
    out[0][0].m = _mm_shuffle_pd(out[0][0].m,out[0][0].m,1);
    out[0][1].m = _mm_shuffle_pd(out[0][1].m,out[0][1].m,1);
    out[0][2].m = _mm_shuffle_pd(out[0][2].m,out[0][2].m,1);

    out[1][0].m = _mm_addsub_pd(psi_[1][0].m, _mm_shuffle_pd(psi_[3][0].m,psi_[3][0].m,1));
    out[1][1].m = _mm_addsub_pd(psi_[1][1].m, _mm_shuffle_pd(psi_[3][1].m,psi_[3][1].m,1));
    out[1][2].m = _mm_addsub_pd(psi_[1][2].m, _mm_shuffle_pd(psi_[3][2].m,psi_[3][2].m,1));

    out[2][0].m = _mm_addsub_pd(psi_[2][0].m, _mm_shuffle_pd(psi_[0][0].m,psi_[0][0].m,1));
    out[2][1].m = _mm_addsub_pd(psi_[2][1].m, _mm_shuffle_pd(psi_[0][1].m,psi_[0][1].m,1));
    out[2][2].m = _mm_addsub_pd(psi_[2][2].m, _mm_shuffle_pd(psi_[0][2].m,psi_[0][2].m,1));

    out[3][0].m = _mm_addsub_pd(_mm_shuffle_pd(psi_[3][0].m,psi_[3][0].m,1), psi_[1][0].m );
    out[3][1].m = _mm_addsub_pd(_mm_shuffle_pd(psi_[3][1].m,psi_[3][1].m,1), psi_[1][1].m );
    out[3][2].m = _mm_addsub_pd(_mm_shuffle_pd(psi_[3][2].m,psi_[3][2].m,1), psi_[1][2].m );
    out[3][0].m = _mm_shuffle_pd(out[3][0].m,out[3][0].m,1);
    out[3][1].m = _mm_shuffle_pd(out[3][1].m,out[3][1].m,1);
    out[3][2].m = _mm_shuffle_pd(out[3][2].m,out[3][2].m,1);
#else
    out[0][0].re = psi_[0][ord][0].re + psi_[2][ord][0].im;    out[0][0].im = psi_[0][ord][0].im - psi_[2][ord][0].re;
    out[0][1].re = psi_[0][ord][1].re + psi_[2][ord][1].im;    out[0][1].im = psi_[0][ord][1].im - psi_[2][ord][1].re;
    out[0][2].re = psi_[0][ord][2].re + psi_[2][ord][2].im;    out[0][2].im = psi_[0][ord][2].im - psi_[2][ord][2].re;
    out[1][0].re = psi_[1][ord][0].re - psi_[3][ord][0].im;    out[1][0].im = psi_[1][ord][0].im + psi_[3][ord][0].re;
    out[1][1].re = psi_[1][ord][1].re - psi_[3][ord][1].im;    out[1][1].im = psi_[1][ord][1].im + psi_[3][ord][1].re;
    out[1][2].re = psi_[1][ord][2].re - psi_[3][ord][2].im;    out[1][2].im = psi_[1][ord][2].im + psi_[3][ord][2].re;
    out[2][0].re = psi_[2][ord][0].re - psi_[0][ord][0].im;    out[2][0].im = psi_[2][ord][0].im + psi_[0][ord][0].re;
    out[2][1].re = psi_[2][ord][1].re - psi_[0][ord][1].im;    out[2][1].im = psi_[2][ord][1].im + psi_[0][ord][1].re;
    out[2][2].re = psi_[2][ord][2].re - psi_[0][ord][2].im;    out[2][2].im = psi_[2][ord][2].im + psi_[0][ord][2].re;
    out[3][0].re = psi_[3][ord][0].re + psi_[1][ord][0].im;    out[3][0].im = psi_[3][ord][0].im - psi_[1][ord][0].re;
    out[3][1].re = psi_[3][ord][1].re + psi_[1][ord][1].im;    out[3][1].im = psi_[3][ord][1].im - psi_[1][ord][1].re;
    out[3][2].re = psi_[3][ord][2].re + psi_[1][ord][2].im;    out[3][2].im = psi_[3][ord][2].im - psi_[1][ord][2].re;
#endif
    break;

  case 3:
#ifdef __INTEL_INTRINSIC__
    out[0][0].m = _mm_add_pd(psi_[0][0].m, psi_[2][0].m);
    out[0][1].m = _mm_add_pd(psi_[0][1].m, psi_[2][1].m);
    out[0][2].m = _mm_add_pd(psi_[0][2].m, psi_[2][2].m);
    out[1][0].m = _mm_add_pd(psi_[1][0].m, psi_[3][0].m);
    out[1][1].m = _mm_add_pd(psi_[1][1].m, psi_[3][1].m);
    out[1][2].m = _mm_add_pd(psi_[1][2].m, psi_[3][2].m);
    out[2][0].m = _mm_add_pd(psi_[2][0].m, psi_[0][0].m);
    out[2][1].m = _mm_add_pd(psi_[2][1].m, psi_[0][1].m);
    out[2][2].m = _mm_add_pd(psi_[2][2].m, psi_[0][2].m);
    out[3][0].m = _mm_add_pd(psi_[3][0].m, psi_[1][0].m);
    out[3][1].m = _mm_add_pd(psi_[3][1].m, psi_[1][1].m);
    out[3][2].m = _mm_add_pd(psi_[3][2].m, psi_[1][2].m);
#else
    out[0][0].re = psi_[0][ord][0].re + psi_[2][ord][0].re;    out[0][0].im = psi_[0][ord][0].im + psi_[2][ord][0].im;
    out[0][1].re = psi_[0][ord][1].re + psi_[2][ord][1].re;    out[0][1].im = psi_[0][ord][1].im + psi_[2][ord][1].im;
    out[0][2].re = psi_[0][ord][2].re + psi_[2][ord][2].re;    out[0][2].im = psi_[0][ord][2].im + psi_[2][ord][2].im;
    out[1][0].re = psi_[1][ord][0].re + psi_[3][ord][0].re;    out[1][0].im = psi_[1][ord][0].im + psi_[3][ord][0].im;
    out[1][1].re = psi_[1][ord][1].re + psi_[3][ord][1].re;    out[1][1].im = psi_[1][ord][1].im + psi_[3][ord][1].im;
    out[1][2].re = psi_[1][ord][2].re + psi_[3][ord][2].re;    out[1][2].im = psi_[1][ord][2].im + psi_[3][ord][2].im;
    out[2][0].re = psi_[2][ord][0].re + psi_[0][ord][0].re;    out[2][0].im = psi_[2][ord][0].im + psi_[0][ord][0].im;
    out[2][1].re = psi_[2][ord][1].re + psi_[0][ord][1].re;    out[2][1].im = psi_[2][ord][1].im + psi_[0][ord][1].im;
    out[2][2].re = psi_[2][ord][2].re + psi_[0][ord][2].re;    out[2][2].im = psi_[2][ord][2].im + psi_[0][ord][2].im;
    out[3][0].re = psi_[3][ord][0].re + psi_[1][ord][0].re;    out[3][0].im = psi_[3][ord][0].im + psi_[1][ord][0].im;
    out[3][1].re = psi_[3][ord][1].re + psi_[1][ord][1].re;    out[3][1].im = psi_[3][ord][1].im + psi_[1][ord][1].im;
    out[3][2].re = psi_[3][ord][2].re + psi_[1][ord][2].re;    out[3][2].im = psi_[3][ord][2].im + psi_[1][ord][2].im;
#endif
    break;
  }

} // uno_p_gmu


template  <int ORD, int DIM>
void BGptSpinColor<ORD,DIM>::uno_m_gmu(SpinColor& out, int mu, int ord){

  switch (mu) {
  case 0:
#ifdef __INTEL_INTRINSIC__
    out[0][0].m = _mm_addsub_pd(psi_[0][0].m, _mm_shuffle_pd(psi_[3][0].m, psi_[3][0].m,1));
    out[0][1].m = _mm_addsub_pd(psi_[0][1].m, _mm_shuffle_pd(psi_[3][1].m, psi_[3][1].m,1));
    out[0][2].m = _mm_addsub_pd(psi_[0][2].m, _mm_shuffle_pd(psi_[3][2].m, psi_[3][2].m,1));
    out[1][0].m = _mm_addsub_pd(psi_[1][0].m, _mm_shuffle_pd(psi_[2][0].m, psi_[2][0].m,1));
    out[1][1].m = _mm_addsub_pd(psi_[1][1].m, _mm_shuffle_pd(psi_[2][1].m, psi_[2][1].m,1));
    out[1][2].m = _mm_addsub_pd(psi_[1][2].m, _mm_shuffle_pd(psi_[2][2].m, psi_[2][2].m,1));

    out[2][0].m = _mm_addsub_pd(_mm_shuffle_pd(psi_[2][0].m,psi_[2][0].m,1), psi_[1][0].m );
    out[2][1].m = _mm_addsub_pd(_mm_shuffle_pd(psi_[2][1].m,psi_[2][1].m,1), psi_[1][1].m );
    out[2][2].m = _mm_addsub_pd(_mm_shuffle_pd(psi_[2][2].m,psi_[2][2].m,1), psi_[1][2].m );
    out[2][0].m = _mm_shuffle_pd(out[2][0].m,out[2][0].m,1);
    out[2][1].m = _mm_shuffle_pd(out[2][1].m,out[2][1].m,1);
    out[2][2].m = _mm_shuffle_pd(out[2][2].m,out[2][2].m,1);

    out[3][0].m = _mm_addsub_pd(_mm_shuffle_pd(psi_[3][0].m,psi_[3][0].m,1), psi_[0][0].m );
    out[3][1].m = _mm_addsub_pd(_mm_shuffle_pd(psi_[3][1].m,psi_[3][1].m,1), psi_[0][1].m );
    out[3][2].m = _mm_addsub_pd(_mm_shuffle_pd(psi_[3][2].m,psi_[3][2].m,1), psi_[0][2].m );
    out[3][0].m = _mm_shuffle_pd(out[3][0].m,out[3][0].m,1);
    out[3][1].m = _mm_shuffle_pd(out[3][1].m,out[3][1].m,1);
    out[3][2].m = _mm_shuffle_pd(out[3][2].m,out[3][2].m,1);
#else
    out[0][0].re = psi_[0][ord][0].re - psi_[3][ord][0].im;      out[0][0].im = psi_[0][ord][0].im + psi_[3][ord][0].re;  
    out[0][1].re = psi_[0][ord][1].re - psi_[3][ord][1].im;      out[0][1].im = psi_[0][ord][1].im + psi_[3][ord][1].re;  
    out[0][2].re = psi_[0][ord][2].re - psi_[3][ord][2].im;      out[0][2].im = psi_[0][ord][2].im + psi_[3][ord][2].re;  
    out[1][0].re = psi_[1][ord][0].re - psi_[2][ord][0].im;      out[1][0].im = psi_[1][ord][0].im + psi_[2][ord][0].re;  
    out[1][1].re = psi_[1][ord][1].re - psi_[2][ord][1].im;      out[1][1].im = psi_[1][ord][1].im + psi_[2][ord][1].re;  
    out[1][2].re = psi_[1][ord][2].re - psi_[2][ord][2].im;      out[1][2].im = psi_[1][ord][2].im + psi_[2][ord][2].re;  
    out[2][0].re = psi_[2][ord][0].re + psi_[1][ord][0].im;      out[2][0].im = psi_[2][ord][0].im - psi_[1][ord][0].re;  
    out[2][1].re = psi_[2][ord][1].re + psi_[1][ord][1].im;      out[2][1].im = psi_[2][ord][1].im - psi_[1][ord][1].re;  
    out[2][2].re = psi_[2][ord][2].re + psi_[1][ord][2].im;      out[2][2].im = psi_[2][ord][2].im - psi_[1][ord][2].re;  
    out[3][0].re = psi_[3][ord][0].re + psi_[0][ord][0].im;      out[3][0].im = psi_[3][ord][0].im - psi_[0][ord][0].re;  
    out[3][1].re = psi_[3][ord][1].re + psi_[0][ord][1].im;      out[3][1].im = psi_[3][ord][1].im - psi_[0][ord][1].re;  
    out[3][2].re = psi_[3][ord][2].re + psi_[0][ord][2].im;      out[3][2].im = psi_[3][ord][2].im - psi_[0][ord][2].re;  
#endif
    break;

  case 1:
#ifdef __INTEL_INTRINSIC__
    out[0][0].m = _mm_add_pd(psi_[0][0].m, psi_[3][0].m);
    out[0][1].m = _mm_add_pd(psi_[0][1].m, psi_[3][1].m);
    out[0][2].m = _mm_add_pd(psi_[0][2].m, psi_[3][2].m);
    out[1][0].m = _mm_sub_pd(psi_[1][0].m, psi_[2][0].m);
    out[1][1].m = _mm_sub_pd(psi_[1][1].m, psi_[2][1].m);
    out[1][2].m = _mm_sub_pd(psi_[1][2].m, psi_[2][2].m);
    out[2][0].m = _mm_sub_pd(psi_[2][0].m, psi_[1][0].m);
    out[2][1].m = _mm_sub_pd(psi_[2][1].m, psi_[1][1].m);
    out[2][2].m = _mm_sub_pd(psi_[2][2].m, psi_[1][2].m);
    out[3][0].m = _mm_add_pd(psi_[3][0].m, psi_[0][0].m);
    out[3][1].m = _mm_add_pd(psi_[3][1].m, psi_[0][1].m);
    out[3][2].m = _mm_add_pd(psi_[3][2].m, psi_[0][2].m);
#else
     out[0][0].re = psi_[0][ord][0].re + psi_[3][ord][0].re;    out[0][0].im = psi_[0][ord][0].im + psi_[3][ord][0].im;
     out[0][1].re = psi_[0][ord][1].re + psi_[3][ord][1].re;    out[0][1].im = psi_[0][ord][1].im + psi_[3][ord][1].im;
     out[0][2].re = psi_[0][ord][2].re + psi_[3][ord][2].re;    out[0][2].im = psi_[0][ord][2].im + psi_[3][ord][2].im;
     out[1][0].re = psi_[1][ord][0].re - psi_[2][ord][0].re;    out[1][0].im = psi_[1][ord][0].im - psi_[2][ord][0].im;
     out[1][1].re = psi_[1][ord][1].re - psi_[2][ord][1].re;    out[1][1].im = psi_[1][ord][1].im - psi_[2][ord][1].im;
     out[1][2].re = psi_[1][ord][2].re - psi_[2][ord][2].re;    out[1][2].im = psi_[1][ord][2].im - psi_[2][ord][2].im;
     out[2][0].re = psi_[2][ord][0].re - psi_[1][ord][0].re;    out[2][0].im = psi_[2][ord][0].im - psi_[1][ord][0].im;
     out[2][1].re = psi_[2][ord][1].re - psi_[1][ord][1].re;    out[2][1].im = psi_[2][ord][1].im - psi_[1][ord][1].im;
     out[2][2].re = psi_[2][ord][2].re - psi_[1][ord][2].re;    out[2][2].im = psi_[2][ord][2].im - psi_[1][ord][2].im;
     out[3][0].re = psi_[3][ord][0].re + psi_[0][ord][0].re;    out[3][0].im = psi_[3][ord][0].im + psi_[0][ord][0].im;
     out[3][1].re = psi_[3][ord][1].re + psi_[0][ord][1].re;    out[3][1].im = psi_[3][ord][1].im + psi_[0][ord][1].im;
     out[3][2].re = psi_[3][ord][2].re + psi_[0][ord][2].re;    out[3][2].im = psi_[3][ord][2].im + psi_[0][ord][2].im;
#endif
     break;

  case 2:
#ifdef __INTEL_INTRINSIC__
    out[0][0].m = _mm_addsub_pd(psi_[0][0].m, _mm_shuffle_pd(psi_[2][0].m, psi_[2][0].m,1));
    out[0][1].m = _mm_addsub_pd(psi_[0][1].m, _mm_shuffle_pd(psi_[2][1].m, psi_[2][1].m,1));
    out[0][2].m = _mm_addsub_pd(psi_[0][2].m, _mm_shuffle_pd(psi_[2][2].m, psi_[2][2].m,1));

    out[1][0].m = _mm_addsub_pd(_mm_shuffle_pd(psi_[1][0].m,psi_[1][0].m,1), psi_[3][0].m );
    out[1][1].m = _mm_addsub_pd(_mm_shuffle_pd(psi_[1][1].m,psi_[1][1].m,1), psi_[3][1].m );
    out[1][2].m = _mm_addsub_pd(_mm_shuffle_pd(psi_[1][2].m,psi_[1][2].m,1), psi_[3][2].m );
    out[1][0].m = _mm_shuffle_pd(out[1][0].m,out[1][0].m,1);
    out[1][1].m = _mm_shuffle_pd(out[1][1].m,out[1][1].m,1);
    out[1][2].m = _mm_shuffle_pd(out[1][2].m,out[1][2].m,1);

    out[2][0].m = _mm_addsub_pd(_mm_shuffle_pd(psi_[2][0].m,psi_[2][0].m,1), psi_[0][0].m );
    out[2][1].m = _mm_addsub_pd(_mm_shuffle_pd(psi_[2][1].m,psi_[2][1].m,1), psi_[0][1].m );
    out[2][2].m = _mm_addsub_pd(_mm_shuffle_pd(psi_[2][2].m,psi_[2][2].m,1), psi_[0][2].m );
    out[2][0].m = _mm_shuffle_pd(out[2][0].m,out[2][0].m,1);
    out[2][1].m = _mm_shuffle_pd(out[2][1].m,out[2][1].m,1);
    out[2][2].m = _mm_shuffle_pd(out[2][2].m,out[2][2].m,1);

    out[3][0].m = _mm_addsub_pd(psi_[3][0].m, _mm_shuffle_pd(psi_[1][0].m, psi_[1][0].m,1));
    out[3][1].m = _mm_addsub_pd(psi_[3][1].m, _mm_shuffle_pd(psi_[1][1].m, psi_[1][1].m,1));
    out[3][2].m = _mm_addsub_pd(psi_[3][2].m, _mm_shuffle_pd(psi_[1][2].m, psi_[1][2].m,1));
#else
    out[0][0].re = psi_[0][ord][0].re - psi_[2][ord][0].im;    out[0][0].im = psi_[0][ord][0].im + psi_[2][ord][0].re;
    out[0][1].re = psi_[0][ord][1].re - psi_[2][ord][1].im;    out[0][1].im = psi_[0][ord][1].im + psi_[2][ord][1].re;
    out[0][2].re = psi_[0][ord][2].re - psi_[2][ord][2].im;    out[0][2].im = psi_[0][ord][2].im + psi_[2][ord][2].re;
    out[1][0].re = psi_[1][ord][0].re + psi_[3][ord][0].im;    out[1][0].im = psi_[1][ord][0].im - psi_[3][ord][0].re;
    out[1][1].re = psi_[1][ord][1].re + psi_[3][ord][1].im;    out[1][1].im = psi_[1][ord][1].im - psi_[3][ord][1].re;
    out[1][2].re = psi_[1][ord][2].re + psi_[3][ord][2].im;    out[1][2].im = psi_[1][ord][2].im - psi_[3][ord][2].re;
    out[2][0].re = psi_[2][ord][0].re + psi_[0][ord][0].im;    out[2][0].im = psi_[2][ord][0].im - psi_[0][ord][0].re;
    out[2][1].re = psi_[2][ord][1].re + psi_[0][ord][1].im;    out[2][1].im = psi_[2][ord][1].im - psi_[0][ord][1].re;
    out[2][2].re = psi_[2][ord][2].re + psi_[0][ord][2].im;    out[2][2].im = psi_[2][ord][2].im - psi_[0][ord][2].re;
    out[3][0].re = psi_[3][ord][0].re - psi_[1][ord][0].im;    out[3][0].im = psi_[3][ord][0].im + psi_[1][ord][0].re;
    out[3][1].re = psi_[3][ord][1].re - psi_[1][ord][1].im;    out[3][1].im = psi_[3][ord][1].im + psi_[1][ord][1].re;
    out[3][2].re = psi_[3][ord][2].re - psi_[1][ord][2].im;    out[3][2].im = psi_[3][ord][2].im + psi_[1][ord][2].re;
#endif

    break;

  case 3:
#ifdef __INTEL_INTRINSIC__
    out[0][0].m = _mm_sub_pd(psi_[0][0].m, psi_[2][0].m);
    out[0][1].m = _mm_sub_pd(psi_[0][1].m, psi_[2][1].m);
    out[0][2].m = _mm_sub_pd(psi_[0][2].m, psi_[2][2].m);
    out[1][0].m = _mm_sub_pd(psi_[1][0].m, psi_[3][0].m);
    out[1][1].m = _mm_sub_pd(psi_[1][1].m, psi_[3][1].m);
    out[1][2].m = _mm_sub_pd(psi_[1][2].m, psi_[3][2].m);
    out[2][0].m = _mm_sub_pd(psi_[2][0].m, psi_[0][0].m);
    out[2][1].m = _mm_sub_pd(psi_[2][1].m, psi_[0][1].m);
    out[2][2].m = _mm_sub_pd(psi_[2][2].m, psi_[0][2].m);
    out[3][0].m = _mm_sub_pd(psi_[3][0].m, psi_[1][0].m);
    out[3][1].m = _mm_sub_pd(psi_[3][1].m, psi_[1][1].m);
    out[3][2].m = _mm_sub_pd(psi_[3][2].m, psi_[1][2].m);
#else
    out[0][0].re = psi_[0][ord][0].re - psi_[2][ord][0].re;    out[0][0].im = psi_[0][ord][0].im - psi_[2][ord][0].im;
    out[0][1].re = psi_[0][ord][1].re - psi_[2][ord][1].re;    out[0][1].im = psi_[0][ord][1].im - psi_[2][ord][1].im;
    out[0][2].re = psi_[0][ord][2].re - psi_[2][ord][2].re;    out[0][2].im = psi_[0][ord][2].im - psi_[2][ord][2].im;
    out[1][0].re = psi_[1][ord][0].re - psi_[3][ord][0].re;    out[1][0].im = psi_[1][ord][0].im - psi_[3][ord][0].im;
    out[1][1].re = psi_[1][ord][1].re - psi_[3][ord][1].re;    out[1][1].im = psi_[1][ord][1].im - psi_[3][ord][1].im;
    out[1][2].re = psi_[1][ord][2].re - psi_[3][ord][2].re;    out[1][2].im = psi_[1][ord][2].im - psi_[3][ord][2].im;
    out[2][0].re = psi_[2][ord][0].re - psi_[0][ord][0].re;    out[2][0].im = psi_[2][ord][0].im - psi_[0][ord][0].im;
    out[2][1].re = psi_[2][ord][1].re - psi_[0][ord][1].re;    out[2][1].im = psi_[2][ord][1].im - psi_[0][ord][1].im;
    out[2][2].re = psi_[2][ord][2].re - psi_[0][ord][2].re;    out[2][2].im = psi_[2][ord][2].im - psi_[0][ord][2].im;
    out[3][0].re = psi_[3][ord][0].re - psi_[1][ord][0].re;    out[3][0].im = psi_[3][ord][0].im - psi_[1][ord][0].im;
    out[3][1].re = psi_[3][ord][1].re - psi_[1][ord][1].re;    out[3][1].im = psi_[3][ord][1].im - psi_[1][ord][1].im;
    out[3][2].re = psi_[3][ord][2].re - psi_[1][ord][2].re;    out[3][2].im = psi_[3][ord][2].im - psi_[1][ord][2].im;
#endif
    break;
  }

}// uno_m_gmu


//template <class B, int AL_ORD, int PT_ORD>
//BGptSU3<B, AL_ORD, PT_ORD> 
//exp( const BGptSU3<B, AL_ORD, PT_ORD>& U ){
//  // Make a copy of U (since exp U = 1 + U + ...)
//  BGptSU3<B, AL_ORD, PT_ORD> result(U);
//  // The pt orders are stored here:
//  pt_q<AL_ORD, PT_ORD> q;
//  // set q = U - V
//  std::copy(result.ptU().begin(), result.ptU().end(),
//            q.ptq().begin());
//  pt_q<AL_ORD, PT_ORD> tmp(q); // store q^i here
//  // calculate \tilde U = exp{q} up to order PT_ORD
//  for (int i = 2; i <= PT_ORD; ++i){
//    tmp *= q; // construct
//    tmp /= i; // q^i / ( i! )
//    result += tmp; // sum it up!
//  }
//  // convert to U^(i) = \tilde U^(i) V 
//  for (int i = 0; i < PT_ORD; ++i)
//    result[i] = result.bgf().ApplyFromRight(result[i]);
//  return result;
//}; 
//
//// needed template specializations...
template void BGptSpinColor<4, 4>::uno_p_gmu(SpinColor&, int, int);
template void BGptSpinColor<4, 4>::uno_m_gmu(SpinColor&, int, int);
//template void BGptSpinColor<6, 4, 4>::uno_p_gmu(SpinColor&, int, int);
//template void BGptSpinColor<6, 4, 4>::uno_m_gmu(SpinColor&, int, int);
//
//template BGptSU3<bgf::AbelianBgf, 6, 2> exp<bgf::AbelianBgf, 6, 2>(const BGptSU3<bgf::AbelianBgf, 6, 2>&);
//template BGptSU3<bgf::AbelianBgf, 6, 4> exp<bgf::AbelianBgf, 6, 4>(const BGptSU3<bgf::AbelianBgf, 6, 4>&);
//template BGptSU3<bgf::AbelianBgf, 10, 8> exp<bgf::AbelianBgf, 10, 8>(const BGptSU3<bgf::AbelianBgf, 10, 8>&);
