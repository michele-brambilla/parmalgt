#ifndef _SU2_ENV_NODE_PT_
#define _SU2_ENV_NODE_PT_


#include <stdlib.h>
#include <stdio.h>
#ifdef __PARALLEL_OMP__
#include <omp.h>
#endif


#include "../MyMath.h"

#include "SU2pt.h"
#include "lattice.h"


#define ADD 1
#define NOADD 0
#define DAG 1
#define NODAG 0

#define DIFF(a, b) (( fabs((a) - (b)) < 0.00000000000001) ? 0 : 1)  
#define DIFF13(a, b) (( fabs((a) - (b)) < 0.0000000000001) ? 0 : 1)  
#define ISDIFF(a, b, p) (( fabs((a) - (b)) < (p)) ? 0 : 1)

#define SGN(a) ((a) >= 0 ? +1 : -1)


#define SWAP(a, b, s) {(s) = (a); (a) = (b); (b) = (s);}


class ptSU2_fld{
 private:
  latt *Z;
 public:
  ptSU2  *W;

  ptSU2_fld(latt* z)  {
    Z = z;
    W = new ptSU2[Z->Size];
  }
  
  ptSU2_fld(FILE*,int); // SU2_fld(&input_file,read_mode)

  ~ptSU2_fld(){
    delete[] W;
  }


  int save(char *filename) {
    FILE *filept;
    if( (filept = fopen(filename, "wb")) == NULL) return 1;
    fwrite(W, sizeof(ptSU2),Z->Size,filept);
    fclose(filept);
    return 0;
  }


  
  int load(char *filename) {
    FILE *filept;
    if( (filept = fopen(filename, "r")) == NULL) return 1;
    fread(W, sizeof(ptSU2),Z->Size,filept);
    fclose(filept);
    return 0;
  }
  

  ptSU2* handle(){ return W; }

  ptSU2 get(int *);

  ptSU2 get(int n){
    return W[Z->L[n][4]];
  }
  
  ptSU2 get(int n,int step,int dir){
    return W[Z->get(n, step, dir)];
  }


  friend int get(ptSU2_fld *, int*);
  friend int get(ptSU2_fld *, int, int, int);
  friend int get(ptSU2_fld *, int);

};
  




inline ptSU2 ptSU2_fld::get(int *n) {
    return W[Z->get(n)];
}



inline int get(ptSU2_fld *W, int *n) {
    return (W->Z)->get(n);

}


inline int get(ptSU2_fld *W,int n,int step,int dir){
  return (W->Z)->get(n,step,dir);
}

inline int get(ptSU2_fld *W,int n){
  return (W->Z)->get(n);
}



// ------------- end class pt SU2_fld  --------------



class ptBoson_fld{
 private:


 public:
  latt   *Z;
  ptBoson  *W;

  ptBoson_fld(latt *z){
    Z = z;
    W = new ptBoson[z->Size];
  }
  
  ptBoson_fld(FILE*,int); // Boson_fld(&input_file,read_mode)


  int save(char *filename) {
    FILE *filept;	
    if( (filept = fopen(filename, "w")) == NULL) return 1;
    if( fwrite(W,sizeof(ptBoson), Z->Size, filept) ) {
      fclose(filept); 
      return 1;
    }
    fclose(filept);
    return 0;
  }

  int load(char *filename) {
    FILE *filept;
    if( (filept = fopen(filename, "r")) == NULL) return 1;
    fread(W,sizeof(ptBoson), Z->Size,filept);
    fclose(filept);
    return 0;
  }


  int save(char *filename, Cplx *placchetta) {
    FILE *filept;	
    if( (filept = fopen(filename, "w")) == NULL) return 1;
    fwrite(W,sizeof(ptBoson), Z->Size, filept);
    fwrite(placchetta, sizeof(Cplx), PTORD+1 , filept);
    fclose(filept);
    return 0;
  }


  int load(char *filename, Cplx *placchetta) {
    FILE *filept;
    if( (filept = fopen(filename, "r") ) == NULL) return 1;
    fread(W,sizeof(ptBoson), Z->Size,filept);
    fread(placchetta, sizeof(Cplx), PTORD+1 , filept);

    fclose(filept);
    return 0;
  }
  

  void operator=(const ptBoson_fld& A) {
    Z = A.Z;
    for(int i = 0; i < Z->Size; i++)
      W[i] = A.W[i];
  }

  ptSU2* handle(){ return (ptSU2*)(W->U); }
  
  ptBoson get(int*);

  ptBoson get(int n, int step, int dir) {return W[Z->get(n, step, dir)];}

  ptSU2 get(int*, int);

  ptSU2 get(int n, int mu){
      return W[Z->get(n)].U[mu];
  }

  ptSU2 get(int n, int step, int dir, int mu){
    return W[Z->get(n, step, dir)].U[mu];
  }
  
  friend int get(ptBoson_fld *W, int n, int mu) {
    return ( 4*((W->Z)->get(n))+mu );
  }

  friend int get(ptBoson_fld *, int *);

  friend int get(ptBoson_fld *W, int n) {
    return ( 4*((W->Z)->get(n)) );
  }

  friend int get(ptBoson_fld *W, int n, int step , int dir, int mu) {
    return (4*((W->Z)->get(n,step,dir))+mu);
  }

  friend int get(ptBoson_fld *W, int n, int step1 , int dir1,
		 int step2, int dir2, int mu) {
    return (4*((W->Z)->get( ((W->Z)->get(n, step1, dir1)), step2, dir2))+mu);
  }


  ptBoson near(int n, int sign, int dir) {
    return W[Z->near(n, sign, dir)-1];
  }    

  ptSU2 staple(int*, int, int); //(x[dim], mu, nu)
  
  inline ptSU2 staple(int, int, int);  //(n, mu, nu)
  inline ptSU2 staple2x1(int, int, int); //(n, mu, nu)
};

  inline ptBoson ptBoson_fld::get(int *n) {
    return W[Z->get(n)];
  }
  
  
  inline ptSU2 ptBoson_fld::get(int *n, int mu) {
    return W[Z->get(n)].U[mu];
  }
  
  
  
inline int get(ptBoson_fld *W, int *n) {
    return W->Z->get(n);
}

ptSU2 ptBoson_fld::staple(int n, int mu, int nu){
  return( W[Z->get(n, 1, mu)].U[nu]*
	  dag(W[Z->L[n][4]].U[nu]*
	      W[Z->get(n, 1, nu)].U[mu]) 
	  +
	  dag(W[Z->get(n, -1, nu)].U[mu]*
	      W[Z->get(n, 1, mu, -1, nu)].U[nu])*
	  W[Z->get(n, -1, nu)].U[nu]);
}

ptSU2 ptBoson_fld::staple2x1(int n, int mu, int nu){

  int curr = Z->L[n][4];
  int a,b,c,d,e,f,g,h;
  ptSU2 tmp;

  /* Contributi alla 2x1. In doppia linea il link in aggiornamento

         d -> 
	 ^    ^
         |    |
    g -> c -> b ->
    ^    ^    ^    ^
    |    |    |    |
    f -> n == a -> e


    h -> n == a -> 
    ^    ^    ^    ^
    |    |    |    |
    g -> c -> b -> f
         ^    ^
	 |    |
	 d -> e

   */

  a = Z->L[curr][5+mu];
  b = Z->L[a][5+nu];
  c = Z->L[curr][5+nu];
  d = Z->L[c][5+nu];
  e = Z->L[a][5+mu];
  f = Z->L[curr][mu  ];
  g = Z->L[f][5+mu];

  tmp = ( W[a].U[nu] * W[b].U[nu] * dag(W[curr].U[nu] * W[c].U[nu] * W[d].U[mu] ) +
          W[a].U[mu] * W[e].U[nu] * dag(W[curr].U[nu] * W[c].U[mu] * W[b].U[mu] ) +
          W[a].U[nu] * dag(W[f].U[nu] * W[g].U[mu] * W[c].U[mu]) * W[f].U[mu] );

  a = Z->L[curr][5+mu];
  b = Z->L[a][nu  ];
  c = Z->L[curr][nu  ];
  d = Z->L[c][nu  ];
  e = Z->L[b][nu  ];
  f = Z->L[b][5+mu];
  g = Z->L[c][mu  ];
  h = Z->L[curr][mu   ];

  tmp += ( dag(W[d].U[mu] * W[e].U[nu] * W[b].U[nu]) * W[d].U[nu] * W[c].U[nu]  +
           W[a].U[mu] * dag(W[c].U[mu] * W[b].U[mu] * W[f].U[nu]) * W[c].U[nu]  +
           dag(W[g].U[mu] * W[c].U[mu] * W[b].U[nu]) * W[g].U[nu] * W[h].U[mu] );

  return tmp;
}

// ********* end class Boson_fld **************

/* class SpinColor_fld { */
    
/* public: */
/*     latt *Z; */
/*     SpinColor *psi; */
  
/*     SpinColor_fld() {} */
/*     void setlatt(latt *z) { */
/* 	Z = z; */
/* 	psi = new SpinColor[Z->Size]; */
/*     } */
    

    
/*     SpinColor_fld(latt *z) { */
/* 	Z = z; */
/* 	psi = new SpinColor[Z->Size]; */
/*     } */

/*     ~SpinColor_fld() { */
/* 	delete[] psi; */
/*     } */

/*     SpinColor get(int n) { */
/* 	return psi[Z->get(n)]; */
/*     } */

/*     SpinColor get(int n, int step, int dir) { */
/* 	return psi[Z->get(n, step, dir)]; */
/*     } */

/*     bool diff(SpinColor_fld *scf) { */
/* 	for (int i = 0; i < Z->Size; i++) */
/* 	    for (int mu = 0; mu < dim; mu++) */
/* 		for (int a = 0; a < NC; a++) */
/* 		    if (DIFF13(psi[Z->get(i)].psi[mu].whr[a].re, scf->psi[Z->get(i)].psi[mu].whr[a].re) || */
/* 			DIFF13(psi[Z->get(i)].psi[mu].whr[a].im, scf->psi[Z->get(i)].psi[mu].whr[a].im)) return (true); */
	
/* 	return (false); */
/*     } */

/*     double norm() { */
/*       double nor = 0.0; */
      
/*       for (int i = 0; i < Z->Size; i++) */
/* 	for (int mu = 0; mu < dim; mu++) */
/* 	  for (int a = 0; a < NC; a++) */
/* 	    nor += ( psi[Z->get(i)].psi[mu].whr[a].re * */
/* 		     psi[Z->get(i)].psi[mu].whr[a].re + */
/* 		     psi[Z->get(i)].psi[mu].whr[a].im * */
/* 		     psi[Z->get(i)].psi[mu].whr[a].im ); */

/*       return (sqrt(nor)); */
/*     } */

/*     void zeros() { */
/* 	for (int i = 0; i < Z->Size; i++) */
/* 	    psi[Z->get(i)].zeros(); */
/*     } */
  
/*     void gauss() { */
/* 	for (int i = 0; i < Z->Size; i++) { */
/* 	    for (int mu = 0; mu < 4; mu++) { */
/* 		for (int c = 0; c < NC; c++) { */
/* 		  psi[i].psi[mu].whr[c].re = generate_gauss(); */
/* 		  psi[i].psi[mu].whr[c].im = generate_gauss(); */
/* 		} */
/* 	    } */
/* 	} */
/*     } */

/*     void base(int *p, int mu, int a) { */
/* 	zeros();//memset(psi, 0, sizeof(SpinColor) * Z->Size); */

/* 	psi[Z->get(p)].psi[mu].whr[a].re = 1.; */
/*     } */

/*     Cplx braket(int *p, int mu, int a) { */
/* 	return (psi[Z->get(p)].psi[mu].whr[a]); */
/*     } */
  
/*     void operator =(SpinColor_fld &scf) { */
/* 	for (int i = 0; i < Z->Size; i++) */
/* 	    for (int mu = 0; mu < dim; mu++) */
/* 		for (int a = 0; a < NC; a++) */
/* 		    psi[i].psi[mu].whr[a] = scf.psi[i].psi[mu].whr[a]; */
/*     } */

/*     void operator +=(SpinColor_fld &scf) { */
/* 	for (int i = 0; i < Z->Size; i++) */
/* 	    for (int mu = 0; mu < dim; mu++) */
/* 		for (int a = 0; a < NC; a++) */
/* 		    psi[i].psi[mu].whr[a] += scf.psi[i].psi[mu].whr[a]; */
/*     } */

/*     void operator -=(SpinColor_fld &scf) { */
/* 	for (int i = 0; i < Z->Size; i++) */
/* 	    psi[i] -= scf.psi[i]; */
/*     } */

/*     void operator *=(double x) { */
/* 	for (int i = 0; i < Z->Size;i++) */
/* 	    psi[i] *= x;	 */
/*     } */
/*     void operator *=(Cplx z) { */
/* 	for (int i = 0; i < Z->Size;i++) */
/* 	    psi[i] *= z;	 */
/*     } */

/*     SpinColor_fld operator *(double x) { */
/* 	SpinColor_fld scf(Z); */
/* 	for (int i = 0; i < Z->Size;i++) */
/* 	    scf.psi[i] = psi[i] * x;	 */
/* 	return scf; */
/*     } */

/*     void prout() { */
/* 	for (int i = 0; i < Z->Size; i++) */
/* 	    for (int mu = 0; mu < dim; mu++) */
/* 		for (int a = 0; a < NC; a++) { */
/* 		    if (DIFF(psi[i].psi[mu].whr[a].re, 0.0) || DIFF(psi[i].psi[mu].whr[a].im, 0.0)) */
/* 			printf("%d %d %d    %+f %+f\n", i, mu, a, psi[i].psi[mu].whr[a].re, psi[i].psi[mu].whr[a].im); */
/* 		} */
/*     } */


/*     SpinColor_fld *pout(dgsp &ds, eigv *ev) { */
/* 	for (int u = 0; u < ds.sz; u++) */
/* 	    psi[Z->get(ev[ds.bs + u].id.i)].psi[ev[ds.bs + u].id.mu].whr[ev[ds.bs + u].id.a] = Cplx(0., 0.); */
	
/* 	return this; */
/*     } */

/*     SpinColor_fld *pin(dgsp &ds, eigv *ev) { */
/* 	for (int u = 0; u < ds.bs; u++) */
/* 	    psi[Z->get(ev[u].id.i)].psi[ev[u].id.mu].whr[ev[u].id.a] = Cplx(0., 0.); */

/* 	for (int u = ds.bs + ds.sz; u < Z->Size*dim*NC; u++) */
/* 	    psi[Z->get(ev[u].id.i)].psi[ev[u].id.mu].whr[ev[u].id.a] = Cplx(0., 0.); */

/* 	return this; */
/*     } */

/*     SpinColor_fld *pin(dgsp &ds, eigv *ev, SpinColor_fld *star) { */
/* 	for (int u = 0; u < ds.bs; u++) */
/* 	    psi[Z->get(ev[u].id.i)].psi[ev[u].id.mu].whr[ev[u].id.a] = Cplx(0., 0.); */

/* 	for (int u = ds.bs + ds.sz; u < Z->Size*dim*NC; u++) */
/* 	    psi[Z->get(ev[u].id.i)].psi[ev[u].id.mu].whr[ev[u].id.a] = Cplx(0., 0.); */

/* 	Cplx c; */
/* 	SpinColor_fld tmp(Z); */
/* 	c = star->scal(this); */
/* 	tmp = *star; */
/* 	tmp *= c; */
/* 	*this -= tmp; */
	
/* 	return this; */
/*     } */

/*     void ket(int i, int mu, int a) { */
/* 	zeros(); */
/* 	psi[Z->get(i)].psi[mu].whr[a].re = 1.; */
/*     } */

/*     void ket(indx id) { */
/* 	zeros(); */
/* 	psi[Z->get(id.i)].psi[id.mu].whr[id.a].re = 1.; */
/*     } */
    
/*     Cplx bra(int i, int mu, int a) { */
/* 	return (psi[i].psi[mu].whr[a]); */
/*     } */

/*     Cplx bra(indx id) { */
/* 	return (psi[id.i].psi[id.mu].whr[id.a]); */
/*     } */

/*     Cplx scal(SpinColor_fld *scf) { */
/* 	Cplx s = Cplx(0, 0); */
/* 	for (int i = 0; i < Z->Size; i++) */
/* 	    for (int mu = 0; mu < dim; mu++) */
/* 		for (int a = 0; a < NC; a++) */
/* 		    s += ~(psi[i].psi[mu].whr[a]) * scf->psi[i].psi[mu].whr[a];  	 */
/* 	return (s); */
/*     } */

/*     SpinColor near(int n, int sign, int dir) { */
/* #ifdef PBC */
/*       return psi[Z->near(n, sign, dir)-1]; */
/* #endif */
/* #ifdef ABC */
/*       int index = Z->near2(n, sign, dir); */
/*       return SGN(index) * psi[abs(index)-1]; */
/* #endif */
/*     } */

/* }; */



/* class ptSpinColor_fld { */
  
/* public: */
/*   latt *Z; */
/*   ptSpinColor *psi; */
  
/*   fftw_plan plan[2]; */
/*   SpinColor_fld *scfld; */

/*   ptSpinColor_fld(latt *z) { */
/*     Z = z; */
    
/*     psi = new ptSpinColor[Z->Size]; */
/*     scfld = new SpinColor_fld(Z); */
    
/*     plan[0] = fftw_plan_many_dft(dim, Z->Sz, NC*dim, */
/* 				 (fftw_complex *) scfld->psi, NULL, NC*dim, 1, */
/* 				 (fftw_complex *) scfld->psi, NULL, NC*dim, 1, */
/* 				 FFTW_FORWARD, FFTW_MEASURE); */
/*     plan[1] = fftw_plan_many_dft(dim, Z->Sz, NC*dim, */
/* 				 (fftw_complex *) scfld->psi, NULL, NC*dim, 1, */
/* 				 (fftw_complex *) scfld->psi, NULL, NC*dim, 1, */
/* 				 FFTW_BACKWARD, FFTW_MEASURE); */
/*   } */
  
  
/*   ~ptSpinColor_fld() { */
/*     delete [] psi; */
/*     delete scfld; */
/*   } */
  
/*   int write(FILE *filept) { */
/*     for (int n = 0; n < Z->Size; n++) */
/*       if (!psi[n].write(filept)) return 0; */
/*     return 1; */
/*   } */
  
/*   int read(FILE *filept) { */
/*     for (int n = 0; n < Z->Size; n++) */
/*       if (!psi[n].read(filept)) return 0;   */
/*     return 1; */
/*   } */
  
/*   int save(char *filename) { */
/*     FILE *filept;	 */
/*     filept = fopen(filename, "w"); */
/*     if( fwrite(psi,sizeof(ptSpinColor), Z->Size,filept) ) { */
/*       fclose(filept); */
/*       return 1; */
/*     } */
/*     fclose(filept); */
/*     return 0; */
/*   } */
  
/*   int load(char *filename) { */
/*     FILE *filept; */
/*     filept = fopen(filename, "r");	 */
/*     if( fread(psi,sizeof(ptSpinColor), Z->Size,filept) ) { */
/*       fclose(filept); */
/*       return 1; */
/*     } */
/*     fclose(filept); */
/*     return 0; */
/*   } */
  
  
/*   void prout(){ */
/*     for(int i = 0; i < Z->Size; i++){ */
/*       printf("i = %d\n",i); */
/*       for(int mu = 0; mu < dim; mu++){ */
/* 	printf("\tmu = %d\n",mu); */
/* 	for(int ord = 0; ord < PTORD+1; ord++){ */
/* 	  printf("\t\tord = %d\n",ord); */
/* 	  psi[i].psi[mu].ptCV[ord].prout(); */
/* 	  printf("\n"); */
/* 	} */
/* 	printf("\n"); */
/*       } */
/*       printf("\n\n"); */
/*     } */
/*   } */
  
/*   void prout(int o){ */
/*     for(int i = 0; i < Z->Size; i++){ */
/*       for(int mu = 0; mu < dim; mu++){ */
/* 	for( int a = 0; a < NC; a++){ */
/* 	  if(fabs(psi[i].psi[mu].ptCV[o].whr[a].re)+ */
/* 	     fabs(psi[i].psi[mu].ptCV[o].whr[a].im) > 1e-10 ){ */
/* 	    printf("%d\t%d\t%d:\t",i,mu,a); */
/* 	    psi[i].psi[mu].ptCV[o].whr[a].prout(); */
/* 	    printf("\n"); */
/* 	  } */
/* 	} */
/*       } */
/*     } */
/*   } */
    
/*   void operator=(ptSpinColor_fld &ptscfld) { */
/*     Z = ptscfld.Z; */
/*     for(int i = 0; i < Z->Size; i++) */
/*       psi[i] = ptscfld.psi[i]; */
/*   } */
  

/*   void scalToPt(int ord) { */
    
/*     for(int i = 0; i < Z->Size; i++){ */
/*       for(int mu = 0; mu < dim; mu++){ */
/* 	psi[i].psi[mu].ptCV[ord] = (scfld->psi[i]).psi[mu]; */
/*       } */
/*     } */

/*   } */
  

/*   ptCVector* handle() { */
/*     return (ptCVector *) (psi->psi); */
/*   } */
  
/*   ptSpinColor get(int *n){ return psi[Z->get(n)]; } */
  
/*   ptSpinColor get(int n, int step, int dir) { */
/*     return psi[Z->get(n, step, dir)]; */
/*   } */
  
/*   ptSpinColor get(int n) { */
/*     return psi[Z->get(n)]; */
/*   } */
  
/*   ptCVector get(int n, int dir) { */
/*     return psi[Z->get(n)].psi[dir]; */
/*   } */


/*   ptSpinColor near(int n, int sign, int dir) { */
/* #ifdef PBC */
/*     return psi[Z->near(n, sign, dir)-1]; */
/* #endif */
/* #ifdef ABC */
/*     int index = Z->near2(n, sign, dir); */
/*     return SGN(index) * psi[abs(index)-1]; */
/* #endif */
/*   } */
  


/*   void fft(int index) { */
    
/* #ifdef ABC */
/*     int coord[4]; */
/*     switch(index){ */
      
/*     case 0: */
/*       for (int i = 0; i < Z->Size; i++){ */
/* 	Z->get(i,coord); */
/* 	scfld->psi[i] *= Cplx(cos(M_PI*coord[0]/Z->Sz[0]), sin(M_PI*coord[0]/Z->Sz[0])); */
/*       } */
/*       fftw_execute(plan[0]); */
/*       break; */
      
/*     case 1: */
/*       fftw_execute(plan[1]); */
/*       for (int i = 0; i < Z->Size; i++){ */
/* 	for (int mu = 0; mu < dim; mu++){ */
/* 	  for (int a = 0; a < NC; a++) { */
/* 	    scfld->psi[i].psi[mu].whr[a].re /= (double) Z->Size; */
/* 	    scfld->psi[i].psi[mu].whr[a].im /= (double) Z->Size; */
/* 	  } */
/* 	} */
	
/* 	Z->get(i,coord); */
/* 	scfld->psi[i] *= Cplx(cos(M_PI*coord[0]/Z->Sz[0]), -sin(M_PI*coord[0]/Z->Sz[0])); */
/*       } */
/*       break; */
/*     } */

/* #elif defined PBC */
/*     fftw_execute(plan[index]); */
/*     if(index==1){ */
/*       for (int i = 0; i < Z->Size; i++){ */
/* 	for (int mu = 0; mu < dim; mu++){ */
/* 	  for (int a = 0; a < NC; a++) { */
/* 	    scfld->psi[i].psi[mu].whr[a] /= (double) Z->Size; */
/* 	  } */
/* 	} */
/*       } */
/*     } */
/* #endif */

/*   } */




/*   void fftT(int index) { */
    
/* #ifdef ABC */
/*     int coord[4]; */
/*     switch(index){ */
      
/*     case 0: */
/*       for (int i = 0; i < Z->Size; i++){ */
/* 	Z->get(i,coord); */
/* 	scfld->psi[i] *= Cplx(cos(M_PI*coord[0]/Z->Sz[0]), sin(M_PI*coord[0]/Z->Sz[0])); */
/*       } */
/*       fftw_execute(plan[1]); */
/*       break; */
      
/*     case 1: */
/*       fftw_execute(plan[0]); */
/*       for (int i = 0; i < Z->Size; i++){ */
/* 	for (int mu = 0; mu < dim; mu++){ */
/* 	  for (int a = 0; a < NC; a++) { */
/* 	    scfld->psi[i].psi[mu].whr[a].re /= (double) Z->Size; */
/* 	    scfld->psi[i].psi[mu].whr[a].im /= (double) Z->Size; */
/* 	  } */
/* 	} */
	
/* 	Z->get(i,coord); */
/* 	scfld->psi[i] *= Cplx(cos(M_PI*coord[0]/Z->Sz[0]), -sin(M_PI*coord[0]/Z->Sz[0])); */
/*       } */
/*       break; */
/*     } */
/* #elif defined PBC */
/*     fftw_execute(plan[(index+1)%2]); */
/*     if(index==1){ */
/*       for (int i = 0; i < Z->Size; i++){ */
/* 	for (int mu = 0; mu < dim; mu++){ */
/* 	  for (int a = 0; a < NC; a++) { */
/* 	    scfld->psi[i].psi[mu].whr[a] /= (double) Z->Size; */
/* 	  } */
/* 	} */
/*       } */
/*     } */
/* #endif */
/*   } */
  
  
/*   void M0inv(){ */
/*     double diag, den; */
/*     Cplx im_(0, -1); */
/*     int i = 0, pp[4]; */
/* #ifdef ABC     */
/*     //    scfld->psi[0] *= Cplx(0.0,0.0); */
/*     for (pp[0] = 0; pp[0] < Z->Sz[0]; pp[0]++) { */
/*       for (pp[1] = 0; pp[1] < Z->Sz[1]; pp[1]++) { */
/* 	for (pp[2] = 0; pp[2] < Z->Sz[2]; pp[2]++) {		 */
/* 	  for (pp[3] = 0; pp[3] < Z->Sz[3]; pp[3]++) { */
/* 	    i = Z->get(pp); */
/* 	    diag = MBARE + 0.5*(Z->p2hat[0][pp[0]] + Z->p2hat[1][pp[1]] +  */
/* 				Z->p2hat[2][pp[2]] + Z->p2hat[3][pp[3]]); */
/* 	    den = diag*diag + (Z->p2bar[0][pp[0]] + Z->p2bar[1][pp[1]] +  */
/* 			       Z->p2bar[2][pp[2]] + Z->p2bar[3][pp[3]]); */
	    
/* 	    scfld->psi[i] = ( ((scfld->psi[i]).gmleft(0) * (Z->pbar[0][pp[0]]) + */
/* 			       (scfld->psi[i]).gmleft(1) * (Z->pbar[1][pp[1]]) + */
/* 			       (scfld->psi[i]).gmleft(2) * (Z->pbar[2][pp[2]]) + */
/* 			       (scfld->psi[i]).gmleft(3) * (Z->pbar[3][pp[3]])) */
/* 			      *im_+ (scfld->psi[i])*diag) / den; */
/* 	  } */
/* 	} */
/*       } */
/*     } */
/* #elif defined PBC */
/*     scfld->psi[Z->get(0)] *= Cplx(0.0,0.0); */
/*     for (pp[0] = 0; pp[0] < Z->Sz[0]; pp[0]++) { */
/*       for (pp[1] = 0; pp[1] < Z->Sz[1]; pp[1]++) { */
/* 	for (pp[2] = 0; pp[2] < Z->Sz[2]; pp[2]++) {		 */
/* 	  for (pp[3] = 0; pp[3] < Z->Sz[3]; pp[3]++) { */
/* 	    if( (pp[0]+pp[1]+pp[2]+pp[3])!= 0){ */
/* 	      i = Z->get(pp); */
/* 	      diag = MBARE + 0.5*(Z->p2hat[0][pp[0]] + Z->p2hat[1][pp[1]] +  */
/* 				  Z->p2hat[2][pp[2]] + Z->p2hat[3][pp[3]]); */
/* 	      den = diag*diag + (Z->p2bar[0][pp[0]] + Z->p2bar[1][pp[1]] +  */
/* 				 Z->p2bar[2][pp[2]] + Z->p2bar[3][pp[3]]); */
	      
/* 	      scfld->psi[i] = ( ((scfld->psi[i]).gmleft(0) * (Z->pbar[0][pp[0]]) + */
/* 				 (scfld->psi[i]).gmleft(1) * (Z->pbar[1][pp[1]]) + */
/* 				 (scfld->psi[i]).gmleft(2) * (Z->pbar[2][pp[2]]) + */
/* 				 (scfld->psi[i]).gmleft(3) * (Z->pbar[3][pp[3]])) */
/* 				*im_+ (scfld->psi[i])*diag) / den; */
/* 	    } */
/* 	  } */
/* 	} */
/*       } */
/*     } */
/* #endif */
/*   } */
  


/*   void M0invT(){ */
  
/*     double diag, den; */
/*     Cplx im_(0, -1); */
/*     int i = 0, pp[4]; */
  
/* #ifdef ABC     */
/*     for (pp[0] = 0; pp[0] < Z->Sz[0]; pp[0]++) { */
/*       for (pp[1] = 0; pp[1] < Z->Sz[1]; pp[1]++) { */
/* 	for (pp[2] = 0; pp[2] < Z->Sz[2]; pp[2]++) {		 */
/* 	  for (pp[3] = 0; pp[3] < Z->Sz[3]; pp[3]++) { */
/* 	    i = Z->get(pp); */
/* 	    diag = MBARE + 0.5*(Z->p2hat[0][pp[0]] + Z->p2hat[1][pp[1]] +  */
/* 				Z->p2hat[2][pp[2]] + Z->p2hat[3][pp[3]]); */
/* 	    den = diag*diag + (Z->p2bar[0][pp[0]] + Z->p2bar[1][pp[1]] +  */
/* 			       Z->p2bar[2][pp[2]] + Z->p2bar[3][pp[3]]); */
	    
/* 	    scfld->psi[i] = ( ((scfld->psi[i]).gmleft(0)*gmT[0] * (Z->pbar[0][pp[0]]) + */
/* 			       (scfld->psi[i]).gmleft(1)*gmT[1] * (Z->pbar[1][pp[1]]) + */
/* 			       (scfld->psi[i]).gmleft(2)*gmT[2] * (Z->pbar[2][pp[2]]) + */
/* 			       (scfld->psi[i]).gmleft(3)*gmT[3] * (Z->pbar[3][pp[3]])) */
/* 			      *im_+ (scfld->psi[i])*diag) / den; */
/* 	  } */
/* 	} */
/*       } */
/*     } */
/* #elif defined PBC */
/*     scfld->psi[Z->get(0)] *= Cplx(0.0,0.0); */
/*     for (pp[0] = 0; pp[0] < Z->Sz[0]; pp[0]++) { */
/*       for (pp[1] = 0; pp[1] < Z->Sz[1]; pp[1]++) { */
/* 	for (pp[2] = 0; pp[2] < Z->Sz[2]; pp[2]++) {		 */
/* 	  for (pp[3] = 0; pp[3] < Z->Sz[3]; pp[3]++) { */
/* 	    if( (pp[0]+pp[1]+pp[2]+pp[3])!= 0){ */
/* 	      i = Z->get(pp); */
/* 	      diag = MBARE + 0.5*(Z->p2hat[0][pp[0]] + Z->p2hat[1][pp[1]] +  */
/* 				  Z->p2hat[2][pp[2]] + Z->p2hat[3][pp[3]]); */
/* 	      den = diag*diag + (Z->p2bar[0][pp[0]] + Z->p2bar[1][pp[1]] +  */
/* 				 Z->p2bar[2][pp[2]] + Z->p2bar[3][pp[3]]); */
	      
/* 	      scfld->psi[i] = ( ((scfld->psi[i]).gmleft(0)*gmT[0] * (Z->pbar[0][pp[0]]) + */
/* 				 (scfld->psi[i]).gmleft(1)*gmT[1] * (Z->pbar[1][pp[1]]) + */
/* 				 (scfld->psi[i]).gmleft(2)*gmT[2] * (Z->pbar[2][pp[2]]) + */
/* 				 (scfld->psi[i]).gmleft(3)*gmT[3] * (Z->pbar[3][pp[3]])) */
/* 				*im_+ (scfld->psi[i])*diag) / den; */
/* 	    } */
/* 	  } */
/* 	} */
/*       } */
/*     } */
/* #endif */
/*   } */





/*   void dirac(int ptord, ptBoson_fld &Umu)  { */
/*     SpinColor S1, S2, S; */
    
/*     for (int i = 0; i < Z->Size; i++) { */
/*       S *= 0;	 */
/* /\*       point_c = Z->get(i); *\/ */

/*       for (int ord = ptord; ord > 0; ord--) { */
	
/* 	for (int nu = 0; nu < dim; nu++) { */
/* 	  S.psi[nu] += (Cplx(mcpt[ord], 0) * get(i).psi[nu].ptCV[ptord-ord]); */
/* /\* 	  S.psi[nu] += (Cplx(mcpt[ord], 0) * psi[point_c].psi[nu].ptCV[ptord-ord]); *\/ */
/* 	} */
	
/* 	for (int mu = 0; mu < dim; mu++) { */
/* /\* 	  point_u = Z->get( i, 1, mu); *\/ */
/* /\* 	  point_d = Z->get( i,-1, mu); *\/ */
	  
/* 	  for (int nu = 0; nu < dim; nu++) { */
/* 	    S1.psi[nu] = Umu.get(i, mu).ptU[ord-1]* near(i, +1, mu).psi[nu].ptCV[ptord-ord]; */
/* 	    S2.psi[nu] = dag(Umu.near(i, -1, mu).U[mu].ptU[ord-1])* near(i, -1, mu).psi[nu].ptCV[ptord-ord]; */
/* /\* 	    S1.psi[nu] = Umu.W[point_c].U[mu].ptU[ord-1] * psi[point_u].psi[nu].ptCV[ptord-ord]; *\/ */
/* /\* 	    S2.psi[nu] = dag(Umu.W[point_d].U[mu]).ptU[ord-1]* psi[point_d].psi[nu].ptCV[ptord-ord]; *\/ */
/* 	  } */

/* 	  S1 -= S1.gmleft(mu); */
/* 	  S2 += S2.gmleft(mu); */
/* 	  S -= (S1 + S2)*.5; */
/* 	} */
	
/*       } */
      
/*       scfld->psi[Z->get(i)] = S; */
      
/*     } */
/*   } */



/* void diracT(int ptord, ptBoson_fld &Umu)  { */
/*     SpinColor S1, S2, S; */
/*     int point_c, point_u, point_d; */
    
/*     for (int i = 0; i < Z->Size; i++) { */
/*       S *= 0;	 */
/*       point_c = Z->get(i); */

/*       for (int ord = ptord; ord > 0; ord--) { */
	
/* 	for (int nu = 0; nu < dim; nu++) { */
/* 	  S.psi[nu] += (Cplx(mcpt[ord], 0) * psi[point_c].psi[nu].ptCV[ptord-ord]); */
/* 	} */
		
/* 	for (int mu = 0; mu < dim; mu++) { */
/* 	  point_u = Z->get( i, 1, mu); */
/* 	  point_d = Z->get( i,-1, mu); */

/* 	  for (int nu = 0; nu < dim; nu++) { */
/* 	    S1.psi[nu] = (~(Umu.W[point_c].U[mu].ptU[ord-1])) * psi[point_u].psi[nu].ptCV[ptord-ord]; */
/* 	    S2.psi[nu] = (tra(Umu.W[point_d].U[mu].ptU[ord-1])) * psi[point_d].psi[nu].ptCV[ptord-ord]; */
/* 	  } */
	  
/* 	  S1 += S1.gmleft(mu) * gmT[mu]; */
/* 	  S2 -= S2.gmleft(mu) * gmT[mu]; */
/* 	  S += (S1 + S2)*.5; */
/* 	} */
	
/*       } */
      
/*       scfld->psi[point_c] = S; */
      
/*     } */
/*   } */






/* /\*   void fillPT(int ptmax, ptBoson_fld &Umu){ *\/ */

/* /\*     fft(0); *\/ */
/* /\*     M0inv(); *\/ */
/* /\*     fft(1); *\/ */
    
/* /\*     for(int i = 0; i < Z->Size; i++){ *\/ */
/* /\*       for(int mu = 0; mu < dim; mu++){ *\/ */
/* /\* 	psi[i].psi[mu].ptCV[0] = scfld->psi[i].psi[mu]; *\/ */
/* /\*       } *\/ */
/* /\*     } *\/ */

/* /\*     for(int ord = 1; ord <= ptmax; ord++){ *\/ */

/* /\*       dirac(ord, Umu); *\/ */
/* /\*       fft(0); *\/ */
/* /\*       M0inv(); *\/ */
/* /\*       fft(1); *\/ */

/* /\*       for(int i = 0; i < Z->Size; i++){ *\/ */
/* /\* 	for(int mu = 0; mu < dim; mu++){ *\/ */
/* /\* 	  psi[i].psi[mu].ptCV[ord] = -(scfld->psi[i].psi[mu]); *\/ */
/* /\* 	} *\/ */
/* /\*       } *\/ */
/* /\*     } *\/ */
/* /\*   } *\/ */



/*   void fillPT(int ptmax, ptBoson_fld &Umu){ */

/*     fft(0); */
/*     M0inv(); */
/*     fft(1); */

/*     int site_c, point_up, point_dn; */
/*     SpinColor Xi1, Xi2; */

/*     // per copiare non mi interessa l'ordine in cui accedo ai siti */
/*     for(int i = 0; i < Z->Size; i++){ */
/*       for(int mu = 0; mu < dim; mu++){ */
/* 	psi[i].psi[mu].ptCV[0] = scfld->psi[i].psi[mu]; */
/*       } */
/*     } */

/*     for(int jord = 1; jord <= ptmax; jord++){ */
      
/*       for(int curr = 0; curr < Z->Size; curr++){ */
	
/* 	site_c = Z->get(curr); */
	      
/* 	// moltiplicazione perturbativa per M */
/* 	for( int kord = 0; kord < jord; kord++) { */
	  
/* 	  // massa critica */
/* 	  for (int mu = 0; mu < dim; mu++ ) {		   */
/* 	    psi[site_c].psi[mu].ptCV[jord] += mcpt[jord-kord]*psi[site_c].psi[mu].ptCV[kord]; */
/* 	  } */
	  
/* 	  for( int mu = 0; mu < dim; mu++){ */
	    
/* 	    point_up = Umu.Z->get(curr, 1, mu); */
/* 	    point_dn = Umu.Z->get(curr,-1, mu); */
	    
/* 	    for( int nu = 0; nu < dim; nu++ ){ */
	      
/* 	      Xi1.psi[nu] = psi[point_dn].psi[nu].ptCV[kord]; */
/* 	      Xi2.psi[nu] = psi[point_up].psi[nu].ptCV[kord]; */
	      
/* 	    } // nu */
	    
/* 	    Xi1 += Xi1.gmleft(mu); */
/* 	    Xi2 -= Xi2.gmleft(mu); */
	    
/* 	    for( int nu = 0; nu < dim; nu++ ){ */
	      
/* 	      psi[site_c].psi[nu].ptCV[jord] -= .5*dag(Umu.W[point_dn].U[mu].ptU[jord-kord-1])*Xi1.psi[nu]; */
/* 	      psi[site_c].psi[nu].ptCV[jord] -= .5*( Umu.W[site_c].U[mu].ptU[jord-kord-1] )*Xi2.psi[nu]; */
/* 	    } //nu */
	    
/* 	  } // mu */
	  
/* 	} //kord */
	
/*       } // curr */
      
      
/*       // metto su scfld per applicare M0^-1, inverto e rimetto in psi */
/*       for(int i = 0; i < Z->Size; i++){ */
/* 	for(int mu = 0; mu < dim; mu++){ */
/* 	  scfld->psi[i].psi[mu] = psi[i].psi[mu].ptCV[jord]; */
/* 	} */
/*       } */
      
/*       fft(0); */
/*       M0inv(); */
/*       fft(1); */
      
/*       for(int i = 0; i < Z->Size; i++){ */
/* 	for(int mu = 0; mu < dim; mu++){ */
/* 	  psi[i].psi[mu].ptCV[jord] = -(scfld->psi[i].psi[mu]); */
/* 	} */
/*       } */
      
      
/*     }// jord */

/*   } */

/* }; */


#endif
