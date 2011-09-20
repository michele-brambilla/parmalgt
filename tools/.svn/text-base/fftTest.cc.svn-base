#include "../QCDenvNODEpt.h"

#ifdef MACOS
#include <mach/mach_time.h>
void mach_absolute_difference(uint64_t end, uint64_t start, struct timespec *tp) {  
  uint64_t difference = end - start;  
  static mach_timebase_info_data_t info = {0,0};  
  
  if (info.denom == 0)  
    mach_timebase_info(&info);  
  
  uint64_t elapsednano = difference * (info.numer / info.denom);  
  
  tp->tv_sec = elapsednano * 1e-9;  
  tp->tv_nsec = elapsednano - (tp->tv_sec * 1e9);  
}  
#endif



int main(int argc, char *argv[])
{


  int *sz = new int[dim];

  if(argc > 1)
    sz[0] = sz[1] = sz[2] = sz[3] = atoi(argv[1]);
  else
    sz[0] = sz[1] = sz[2] = sz[3] = 8;
  
  latt LL(sz);
  LL.p_init();

  ptSpinColor_fld Psi(&LL);


  {
    MyRand Rand;
    Rand.init(1234);

#ifndef MACOS
  struct timespec t1,t2;
  clock_gettime(CLOCK_REALTIME,&t1);
#else
  uint64_t start,end;  
  struct timespec tp;  
  start = mach_absolute_time();  
#endif

  Psi.scfld->gauss(Rand);
    
  //  Rand.generate_gauss(&Psi.scfld->psi[0].psi[0].whr[0].re,LL.Size*12*2);

#ifndef MACOS
  clock_gettime(CLOCK_REALTIME,&t2);
  printf("t[ms] = %g\n", 1e3*(t2.tv_sec-t1.tv_sec)+1e-6*(t2.tv_nsec-t1.tv_nsec) );
#else
  end = mach_absolute_time();  
  mach_absolute_difference(end, start, &tp);
  printf("Random source generation: t[ms] = %g\n", 1e3*tp.tv_sec+1e-6*tp.tv_nsec);
#endif


  }


#ifndef MACOS
  struct timespec t1,t2;
  clock_gettime(CLOCK_REALTIME,&t1);
#else
  uint64_t start,end;  
  struct timespec tp;  
  start = mach_absolute_time();  
#endif

  Psi.fft(0);
  Psi.fft(1);


#ifndef MACOS
  clock_gettime(CLOCK_REALTIME,&t2);
  printf("2 x FFT: t[ms] = %g\n", 1e3*(t2.tv_sec-t1.tv_sec)+1e-6*(t2.tv_nsec-t1.tv_nsec) );
#else
  end = mach_absolute_time();  
  mach_absolute_difference(end, start, &tp);
  printf("2 x FFT: t[ms] = %g\n", 1e3*tp.tv_sec+1e-6*tp.tv_nsec);
#endif



#ifndef MACOS
  clock_gettime(CLOCK_REALTIME,&t1);
#else
  start = mach_absolute_time();  
#endif

  Psi.fft(0);
  Psi.M0inv();
  Psi.fft(1);


#ifndef MACOS

  clock_gettime(CLOCK_REALTIME,&t2);
  printf("2 x FFT + M0inv: t[ms] = %g\n", 1e3*(t2.tv_sec-t1.tv_sec)+1e-6*(t2.tv_nsec-t1.tv_nsec) );

#else

  end = mach_absolute_time();  
  mach_absolute_difference(end, start, &tp);
  printf("2 x FFT + M0inv: t[ms] = %g\n", 1e3*tp.tv_sec+1e-6*tp.tv_nsec);

#endif


  return 0;

}
