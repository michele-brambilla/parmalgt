#include <cmath>
#include <limits.h>
#include <cstdlib>

#include "dSFMT/dSFMT.h"

#define __GEN_RAND55__
//#define __GEN_RAND__
//#define __GEN_DSFMT__


#define MY_2pi 6.28318530717959 

#ifdef __GEN_RAND55__
#ifdef ULONG_MAX
const double Norm_Factor = 1./(double)ULONG_MAX;
#else
const double Norm_Factor = 2.328306437080797e-10;
#endif
#elif defined __GEN_RAND__
const double Norm_Factor = 1./(double)RAND_MAX;
#endif


class MyRand {
private:
// per il generatore piatto "di base"

  /* FILE* fp; */

#ifdef __GEN_RAND55__
  unsigned long ring[55];
  int next[55];
  int n1, n2;
#elif defined __GEN_DSFMT__
  dsfmt_t mydsfmt;
#endif

  double value;

// per il generatore gaussiano
  bool gen;
  double r,t;

// per il generatore gaussiano di media e larghezza variabile
  int gss_flag;
  double gss_value;
  double gss_sigma;
  double gss_x0;

// funzioni per il generatore piatto 

  double Ranf()  {
#ifdef __GEN_RAND55__
    ring[n1=next[n1]] += ring[n2=next[n2]]; 
    value = (double)ring[n1] * Norm_Factor; 
#elif defined __GEN_RAND__
    value = rand()*Norm_Factor;
#elif defined __GEN_DSFMT__
    value = dsfmt_genrand_close_open(&mydsfmt);
#endif
    return value;
  }




  double* Ranf(double *v, int n)  {

    for( int i = 0; i < n; ++i)
      {
#ifdef __GEN_RAND55__
	ring[n1=next[n1]] += ring[n2=next[n2]]; 
	v[i] = (double)ring[n1] * Norm_Factor; 
#elif defined __GEN_RAND__
	v[i] = rand()*Norm_Factor;
#endif
      }
    return v;
  }

  


  // funzioni per il generatore gaussiano 
  double get_gss(double sgm, double xM) {
    double r1, r2, theta, rho, Xout, Xgss;

    r1 = Ranf();
    r2 = Ranf();

    rho = sqrt(-2.0 * log(r1));
    rho = sgm * rho;

    theta = MY_2pi * r2;

    Xout = rho * cos(theta);
    Xgss = rho * sin(theta);

    Xout = Xout + xM;
    Xgss = Xgss + xM;

    gss_value = Xgss;
    gss_sigma = sgm;
    gss_x0 = xM;
    gss_flag = 1;

    return Xout;
  }


  double get_from_gss(double sgm, double xM){
    gss_flag = 0; 
    if ( (sgm==gss_sigma)&&(xM==gss_x0) ) 
      return gss_value; 
    else
      return ((gss_value-gss_x0)*(sgm/gss_sigma)+xM);
  }
  

public:

  MyRand() { gen = 0; };
  MyRand(long);

  ~MyRand() {};

  void init(long);

  void init(long, int);

  double Rand() {
    Ranf(); 
    return value;
  }

  double* Rand(double* v, int n) {
    return Ranf(v, n); 
  }

  double randN(double sgm, double xM) {
    double tmp; 
    tmp = (gss_flag==0) ? get_gss(sgm,xM) : get_from_gss(sgm,xM); 
    return tmp;
  }

  double generate_gauss();
  /* double generate_gauss(int); */
  double* generate_gauss(double* v, int n);
  
};

