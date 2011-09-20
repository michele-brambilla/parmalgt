#include "MyRand.h"
#include <stdio.h>
#include <string.h>

long Seed = 161803398L;

MyRand::MyRand(long idum)
{
  gen = 0;
#ifdef __GEN_RAND55__
    int i,j;
    long tmp,aux;

    for(i=0;i<54;++i) {
      next[i] = i+1;
    }
    next[54]= 0;
    tmp = Seed + idum;
    ring[54] = tmp;
    aux = 1;
    j = 20;
    for(i=0;i<54;++i,j+=21){
	j %= 55;
	ring[j] = aux;
	aux += tmp;
	tmp = ring[j];
    }

    n1 = 0;
    n2 = 31;

    for(j=0;j<55*4;++j) {
      Ranf();
    }

    gss_flag = 0;
    gss_value = 0.0;
    gss_sigma = 0.0;
    gss_x0 = 0.0;
#elif defined __GEN_RAND__
    srand(idum);
#elif defined __GEN_DSFMT__
    dsfmt_init_gen_rand ( &mydsfmt, idum );
#endif
}


void MyRand::init(long idum)
{
#ifdef __GEN_RAND55__
    int i,j;
    long tmp,aux;

    for(i=0;i<54;++i) {
      next[i] = i+1;
    }
    next[54]= 0;
    tmp = Seed + idum;
    ring[54] = tmp;
    aux = 1;
    j = 20;
    for(i=0;i<54;++i,j+=21){
	j %= 55;
	ring[j] = aux;
	aux += tmp;
	tmp = ring[j];
    }

    n1 = 0;
    n2 = 31;

    for(j=0;j<55*4;++j) {
      Ranf();
    }

    gss_flag = 0;
    gss_value = 0.0;
    gss_sigma = 0.0;
    gss_x0 = 0.0;
#elif defined __GEN_RAND__   
    srand(idum);
#elif defined __GEN_DSFMT__
    dsfmt_init_gen_rand ( &mydsfmt, idum );
#endif
}




double MyRand::generate_gauss(){
  double x;

  if(gen == 0){
    t = MY_2pi*Ranf();
    r = sqrt( -log((1. - Ranf())) );
    x = r*cos(t);
    gen = 1;
  }
  else{
    x = r*sin(t);
    gen = 0;
  }
  
  return(x);
}






double* MyRand::generate_gauss(double* v, int n){
  register double r,t;

  Ranf(v,n);

  for( int i = 0; i < n; i+=2)
    {
      t = MY_2pi*v[i];
      r = sqrt( -log((1. - v[i+1])) );
      v[i]   = r*cos(t);
      v[i+1] = r*sin(t);
  }
  
  return v;
}

