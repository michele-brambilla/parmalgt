#include<iostream>
#include<fstream>

#include<ctime>

#include "input.h"
extern nspt_params_t nspt_pars;

using namespace std;

#ifndef MACOS
void mach_absolute_difference(uint64_t, uint64_t, struct timespec*);
#endif

class PRlgtTime{
 private:
  ofstream of_timing;

  struct timespec t_g,  t_f,  t_zm,  t_gf;

#ifndef MACOS
  struct timespec t0_g, t0_f, t0_zm, t0_gf;
  struct timespec t1_g, t1_f, t1_zm, t1_gf;
#else
  uint64_t start_g, end_g, start_f, end_f, start_zm, end_zm, start_gf, end_gf;
  struct timespec ts_g, ts_f, ts_zm, ts_gf;
#endif

 public:
  
  PRlgtTime() {};
  
  ~PRlgtTime()
    {
      of_timing.close();
    } 

  void init();

  /* tic: start time */
  /* toc: end time   */
  /* for the different parts of the action: */
  /* g -> gluonic action */
  /* f -> fermionic action */
  /* zm -> zero momentum subtraction */
  /* gf -> gauge fixing */

  void tic_g();
  void toc_g();
  void tic_f();
  void toc_f();
  void tic_zm();
  void toc_zm();
  void tic_gf();
  void toc_gf();

  /* Computes the timediffrence and accumulates */
  void reduce();

  /* Print out the time */
  void out();

  /* Reset counters */
  void reset();
};
