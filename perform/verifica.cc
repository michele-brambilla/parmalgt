#undef RANF

#include <iostream>
#include<stdlib.h>
#include<fstream>

#include<cmath>

#include <time.h>
#include<sys/time.h>
#include<ctime>

#include <emmintrin.h> 

#include<QCDpt.h>


const int START = 5; // log10 of start nop
const int STOP  = 9; // log10 of stop  nop


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


typedef struct{
  int i;
  double flops_gettimeofday;
  double flops_clock;
  double flops_clock_gettime;
  
  double tempo_gettimeofday;
  double tempo_clock;
  double tempo_clock_gettime;
} data_t;


data_t tempi;


clock_t start_c, end_c;
timeval tg1, tg2, tg;
#ifndef MACOS
struct timespec ts1,ts2;
#else
  uint64_t start_g,end_g;  
  struct timespec ts; 
#endif

using namespace std;


void stampa_tempi(long nop)
{

  tempi.tempo_gettimeofday  = 1e3*(tg2.tv_sec - tg1.tv_sec) + 1e-3*(tg2.tv_usec - tg1.tv_usec);
  tempi.tempo_clock         = 1e3*(end_c - start_c)/(double)CLOCKS_PER_SEC;
#ifndef MACOS
  tempi.tempo_clock_gettime = 1e3*(ts2.tv_sec - ts1.tv_sec) + 1e-6*(ts2.tv_nsec - ts1.tv_nsec);
#else
  mach_absolute_difference(end_g, start_g, &ts);
  tempi.tempo_clock_gettime = 1e3*ts.tv_sec + 1e-6*ts.tv_nsec;
#endif

  cout << tempi.tempo_gettimeofday      << "\t"
       << nop/tempi.tempo_gettimeofday  << "\t\t";

  cout << tempi.tempo_clock             << "\t"
       << nop/tempi.tempo_clock         << "\t\t";

  cout << tempi.tempo_clock_gettime     << "\t"
       << nop/tempi.tempo_clock_gettime << endl;

}




int main(int argc, char *argv[]) {
  
  cout << "gettimeofday\t\tclock\t\t\tclock_gettime" << endl;

  cout << "ms\tMFps\t\t"
       << "ms\tMFps\t\t"
       << "ms\tMFps" << endl << endl;

  Cplx a(2.345456,1.3453546),b(2.345457,1.3453545),c;
  
  
  
  int x = 5;
  
  gettimeofday(&tg1, NULL);
  start_c = clock();
#ifndef MACOS
  clock_gettime(CLOCK_REALTIME,&ts1);
#else
  start_g = mach_absolute_time();  
#endif
  
  for( uint64_t i = 0; i < 1e5; i++ )
    {
      c += a;
      c -= b;
    }
  
#ifndef MACOS
  clock_gettime(CLOCK_REALTIME,&ts2);
#else
  end_g = mach_absolute_time();  
#endif
  
  gettimeofday(&tg2, NULL);
  end_c = clock();
  
  // Reduce by a factor 1e-3 nop to get correct normalization
  stampa_tempi(4*1e2);      
      

  c.prout();
  cout << endl;





  return 0;
}
