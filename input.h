#ifndef _INPUT_H_
#define _INPUT_H_

#include<iostream>
#include<fstream>

#include"QCDenvNODEpt.h"

#ifdef MACOS
#include <mach/mach_time.h>
#endif


typedef struct{
  int Sweep;
  int Beat;
  int Kill;
  int Save;
  int Init;
  char *plaqn;
  char *confn;
  char *damon;
  char *normn;
  char *trun;
  char *logn;
  char *name2x2;
} nspt_params_t;


typedef struct{
  int iVol;
  int *sz;
  double rVol;
  double tau_f;
  double tau_g;
  double stau;
  double alpha;
  double c0;
  double c1;
} act_params_t;

typedef struct{
  int *xi;
  int *xf;
} thread_params_t;

#endif
