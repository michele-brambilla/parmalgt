#ifndef _REN_LIB_H_
#define _REN_LIB_H_

#include "../QCDenvNODEpt.h"
#include <time.h>


#ifdef _AFFINITY
#include <sys/sysinfo.h>
#include <sched.h>
#include <procps/readproc.h>
#endif

#ifdef __RENORM_OMP__
#include<omp.h>
#endif

#include<fstream>


#define GF_PREC 0.0025


typedef struct{
  int iVol;
  int ptord;
  double alpha;
  int *sz;
  char *conf;
  char *out; 
  char *prfile;
#ifdef __TRANSP_PROPAG__
  char *prfileT;
#endif
  char *mom;
  std::ofstream momout;
} renorm_params_t;


int carica_campo(ptGluon_fld&, Cplx*);

int gauge_fixing(ptGluon_fld&, int&);

void fft(ptSpinColor_fld** , int, int);

void XiBuild1(ptGluon_fld&, void**, int, int);



#endif
