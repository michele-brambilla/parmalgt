#include "MyTime.h"

using namespace std;


#ifdef MACOS
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


void PRlgtTime::init()
{
  char* name = new char[110];
  sprintf(name, "time_%s", nspt_pars.logn);

  of_timing.open(name);
#ifdef __TIMING__
  of_timing << "Gauge\tFermion\tGaugeFix\tZeroMom" << endl << endl; 
#endif

  
  delete [] name;
}


void PRlgtTime::reset()
{
  t_g.tv_sec   = 0;
  t_g.tv_nsec  = 0;
  t_f.tv_sec   = 0;
  t_f.tv_nsec  = 0;
  t_gf.tv_sec  = 0;
  t_gf.tv_nsec = 0;
  t_zm.tv_sec  = 0;
  t_zm.tv_nsec = 0;
}


void PRlgtTime::tic_g()
{
#ifndef MACOS
  clock_gettime(CLOCK_REALTIME,&t0_g);
#else
  start_g = mach_absolute_time();
#endif
}

void PRlgtTime::toc_g()
{
#ifndef MACOS
  clock_gettime(CLOCK_REALTIME,&t1_g);
#else
  end_g   = mach_absolute_time();
#endif
}


void PRlgtTime::tic_f()
{
#ifndef MACOS
  clock_gettime(CLOCK_REALTIME,&t0_f);
#else
  start_f = mach_absolute_time();
#endif
}

void PRlgtTime::toc_f()
{
#ifndef MACOS
  clock_gettime(CLOCK_REALTIME,&t1_f);
#else
  end_f   = mach_absolute_time();
#endif
}



void PRlgtTime::tic_zm()
{
#ifndef MACOS
  clock_gettime(CLOCK_REALTIME,&t0_zm);
#else
  start_zm = mach_absolute_time();
#endif
}

void PRlgtTime::toc_zm()
{
#ifndef MACOS
  clock_gettime(CLOCK_REALTIME,&t1_zm);
#else
  end_zm = mach_absolute_time();
#endif
}


void PRlgtTime::tic_gf()
{
#ifndef MACOS
  clock_gettime(CLOCK_REALTIME,&t0_gf);
#else
  start_gf = mach_absolute_time();
#endif
}

void PRlgtTime::toc_gf()
{
#ifndef MACOS
  clock_gettime(CLOCK_REALTIME,&t1_gf);
#else
  end_gf = mach_absolute_time();
#endif
}



void PRlgtTime::reduce()
{
#ifndef MACOS
    t_g.tv_sec   += (t1_g.tv_sec  - t0_g.tv_sec);
    t_g.tv_nsec  += (t1_g.tv_nsec - t0_g.tv_nsec);

    t_f.tv_sec   += (t1_f.tv_sec  - t0_f.tv_sec);
    t_f.tv_nsec  += (t1_f.tv_nsec - t0_f.tv_nsec);

    t_gf.tv_sec  += (t1_gf.tv_sec  - t0_gf.tv_sec);
    t_gf.tv_nsec += (t1_gf.tv_nsec - t0_gf.tv_nsec);

    t_zm.tv_sec  += (t1_zm.tv_sec  - t0_zm.tv_sec);
    t_zm.tv_nsec += (t1_zm.tv_nsec - t0_zm.tv_nsec);
#else
    mach_absolute_difference(end_g,  start_g,  &ts_g );
    t_g.tv_sec  += ts_g.tv_sec;
    t_g.tv_nsec += ts_g.tv_nsec;

    mach_absolute_difference(end_f,  start_f,  &ts_f );
    t_f.tv_sec  += ts_f.tv_sec;
    t_f.tv_nsec += ts_f.tv_nsec;

    mach_absolute_difference(end_gf, start_gf, &ts_gf);
    t_gf.tv_sec  += ts_gf.tv_sec;
    t_gf.tv_nsec += ts_gf.tv_nsec;

    mach_absolute_difference(end_zm, start_zm, &ts_zm);
    t_zm.tv_sec  += ts_zm.tv_sec;
    t_zm.tv_nsec += ts_zm.tv_nsec;
#endif
}


void PRlgtTime::out()
{
  of_timing << 1e3*t_g.tv_sec  + 1e-6*t_g.tv_nsec  << "\t"
	    << 1e3*t_f.tv_sec  + 1e-6*t_f.tv_nsec  << "\t"
	    << 1e3*t_gf.tv_sec + 1e-6*t_gf.tv_nsec << "\t"
	    << 1e3*t_zm.tv_sec + 1e-6*t_zm.tv_nsec << "\t"
	    << endl;
}
