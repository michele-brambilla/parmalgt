#ifndef SF
#include "QCDpt.h"
#else
#include "newQCDpt.h"
const int ORD = 6;
typedef BGptSU3<bgf::AbelianBgf, ORD> ptSU3;
typedef ptt::PtMatrix<ORD> ptsu3;
typedef BGptCVector<ORD> ptCVector;
typedef BGptGluon<bgf::AbelianBgf, ORD, 4> ptGluon;
typedef BGptSpinColor<ORD, 4> ptSpinColor;
#endif
#include "lattice.h"
#include<iostream>
#include<fstream>


#define DIFF(a, b) (( fabs((a) - (b)) < 1e-15) ? 0 : 1)  
#define DIFF13(a, b) (( fabs((a) - (b)) < 1e-13) ? 0 : 1)  
#define ISDIFF(a, b, p) (( fabs((a) - (b)) < (p)) ? 0 : 1)

using namespace std;


typedef struct {
     int i, mu, a;
} indx;

typedef struct {
    indx id;
    int mo;
    double pt[ORD+1];
    int cf[ORD+1];
} eigv;


class dgsp;


class SU3_fld;
class Gluon_fld;
class SpinColor_fld;
  
  
class SU3_fld{

 public:

  latt *Z;
  SU3  *W;

  SU3_fld(latt* z)  {
    Z = z;
    W = new SU3[Z->Size];
  }

  ~SU3_fld(){
    delete [] W;
  }
  
  SU3_fld(FILE*,int); // SU3_fld(&input_file,read_mode)

  SU3* handle(){ return W; }

  SU3 get(int *);

  const SU3& operator[](int i) const { return W[i]; };
  SU3& operator[](int i) { return W[i]; };
  
  SU3 get(int n,int step,int dir){
    return W[Z->get(n, step, dir)];
  }


  friend int get(SU3_fld *, int*);
  friend int get(SU3_fld *, int, int, int);
};  
  




/*   inline SU3 SU3_fld::get(int *n) { */
/* #if (dim == 2) */
/*     return W[ Z->L[*n][*(n+1)] ]; */
/* #endif */
/* #if (dim == 3) */
/*     return W[ Z->L[*n][*(n+1)][*(n+2)] ]; */
/* #endif */
/* #if (dim == 4) */
/*     return W[ Z->L[*n][*(n+1)][*(n+2)][*(n+3)] ]; */
/* #endif     */
/*   } */



/*   inline int get(SU3_fld *W, int *n) { */
/* #if (dim == 2) */
/*     return W->Z->L[*n][*(n+1)]; */
/* #endif */
/* #if (dim == 3) */
/*     return W->Z->L[*n][*(n+1)][*(n+2)]; */
/* #endif */
/* #if (dim == 4) */
/*     return W->Z->L[*n][*(n+1)][*(n+2)][*(n+3)]; */
/* #endif     */
/*   } */

  
inline int get(SU3_fld *W,int n,int step,int dir){
  return W->Z->get(n,step,dir);
}

  
class Gluon_fld{
 private:
  latt   *Z;

 public:
  Gluon  *W;
  Gluon_fld(latt *z){
    Z = z;
    W = new Gluon[z->Size];
  }
  
  ~Gluon_fld(){
    delete [] W;
  }
  
  Gluon_fld(FILE*,int); // Gluon_fld(&input_file,read_mode)
  

  void operator=(const Gluon_fld& A) {
    Z = A.Z;
    for(int i = 0; i < Z->Size; i++)
      W[i] = A.W[i];
  }
  
  



/*   SU3* handle(){ return (SU3*)W->U; } */
/*   Gluon* _handle(){ return W; } */
  
/*   Gluon get(int*); */
  Gluon get(int n, int step, int dir) {return W[Z->get(n, step, dir)];}
  
/*   SU3 get(int*, int); */
/*   SU3 get(int n, int step, int dir, int mu){ */
/*     SU3 *ptr; */
/*     ptr = (SU3*)(W+Z->get(n, step, dir)); */
/*     ptr += mu; */
/*     return *ptr; */
/*   } */
  
  
/*   // for gluon */
/*   friend int get(Gluon_fld *W, int n, int step, int dir){ */
/*     return W->Z->get(n,step,dir); */
/*   } */
/*   // for SU3* */
/*   friend int get(Gluon_fld *W, int n, int step , int dir, int mu) { */
/*     return 4*(W->Z->get(n,step,dir))+mu; */
/*   } */


/*   friend int get(Gluon_fld *W, int n, int step1 , int dir1,  */
/* 		 int step2, int dir2, int mu) { */
/*     return 4*(W->Z->get(n,step1,dir1, step2, dir2))+mu; */
/*   } */

  
/*   friend int get(Gluon_fld *, int *); */

  
/*   SU3 staple(int*, int, int); //(x[dim], mu, nu) */
  
/*   SU3 staple(int, int, int);  //(n, mu, nu) */
  
/* }; */



/* inline Gluon Gluon_fld::get(int *n) { */
/* #if (dim == 2) */
/*   return W[ Z->L[*n][*(n+1)] ]; */
/* #endif */
/* #if (dim == 3) */
/*   return W[ Z->L[*n][*(n+1)][*(n+2)] ]; */
/* #endif */
/* #if (dim == 4) */
/*   return W[ Z->L[*n][*(n+1)][*(n+2)][*(n+3)] ]; */
/* #endif     */
/* } */


/* inline SU3 Gluon_fld::get(int *n, int mu) { */
/*   Gluon *ptg; */
/*   SU3* ptr; */
/* #if (dim==2) */
/*   ptg = (W + Z->L[*n][*(n+1)]); */
/*   ptr = ((SU3*)ptg) + mu; */
/* #endif */
/* #if (dim==3) */
/*   ptg = (W + Z->L[*n][*(n+1)][*(n+2)]); */
/*   ptr = ((SU3*)ptg) + mu; */
/* #endif */
/* #if (dim==4) */
/*   ptg = (W + Z->L[*n][*(n+1)][*(n+2)][*(n+3)]); */
/*   ptr = ((SU3*)ptg) + mu; */
/* #endif     */
/*   return *ptr; */
/*   } */



/*   inline int get(Gluon_fld *W, int *n) { */
/* #if (dim == 2) */
/*     return W->Z->L[*n][*(n+1)]; */
/* #endif */
/* #if (dim == 3) */
/*     return W->Z->L[*n][*(n+1)][*(n+2)]; */
/* #endif */
/* #if (dim == 4) */
/*     return W->Z->L[*n][*(n+1)][*(n+2)][*(n+3)]; */
/* #endif     */
/*   } */







  SU3 staple(int n, int mu, int nu){

  /*
    b ->
    ^    ^
    |    |
    n == a
    ^    ^
    |    |
    d -> c

   */
  
    return( W[Z->get(n, 1, mu)].U[nu]*
	    dag(W[Z->L[n][4]].U[nu]*
		W[Z->get(n, 1, nu)].U[mu]) 
	    +
	    dag(W[Z->get(n, -1, nu)].U[mu]*
		W[Z->get(n, 1, mu, -1, nu)].U[nu])*
	    W[Z->get(n, -1, nu)].U[nu]);
  }
  
};






class dgsp{
 public:
  double ev;
  int bs, sz, ns;
  dgsp *ds;
  
  dgsp(){};

  int n_subs(eigv*, int);
  
  void build_subs(eigv*, int);

};



class SpinColor_fld {
    
public:
    latt *Z;
    SpinColor *psi;
  
    SpinColor_fld() {}
    void setlatt(latt *z) {
	Z = z;
	psi = new SpinColor[Z->Size];
    }
    
    const SpinColor& operator[](int i) const { return psi[i]; }
    SpinColor& operator[](int i) { return psi[i]; }
    
    
    SpinColor_fld(latt *z) {
	Z = z;
	psi = new SpinColor[Z->Size];
    }

    ~SpinColor_fld() {
	delete[] psi;
    }

    SpinColor get(int n) {
	return psi[Z->get(n)];
    }

    SpinColor get(int n, int step, int dir) {
	return psi[Z->get(n, step, dir)];
    }

    bool diff(SpinColor_fld *scf) {
	for (int i = 0; i < Z->Size; i++)
	    for (int mu = 0; mu < dim; mu++)
		for (int a = 0; a < NC; a++)
		    if (DIFF13(psi[Z->get(i)].psi[mu].whr[a].re, scf->psi[Z->get(i)].psi[mu].whr[a].re) ||
			DIFF13(psi[Z->get(i)].psi[mu].whr[a].im, scf->psi[Z->get(i)].psi[mu].whr[a].im)) return (true);
	
	return (false);
    }

    double norm() {
      double nor = 0.0;
      
      for (int i = 0; i < Z->Size; i++)
	for (int mu = 0; mu < dim; mu++)
	  for (int a = 0; a < NC; a++)
	    nor += ( psi[Z->get(i)].psi[mu].whr[a].re *
		     psi[Z->get(i)].psi[mu].whr[a].re +
		     psi[Z->get(i)].psi[mu].whr[a].im *
		     psi[Z->get(i)].psi[mu].whr[a].im );

      return (sqrt(nor));
    }

    void zeros() {
      bzero(psi, Z->Size*sizeof(SpinColor));
      /* for (int i = 0; i < Z->Size; i++) */
      /* 	psi[Z->get(i)].zeros(); */
    }

    void gauss(MyRand& Rand) {
	for (int i = 0; i < Z->Size; i++) {
	    for (int mu = 0; mu < dim; mu++) {
		for (int c = 0; c < NC; c++) {
		  psi[i].psi[mu].whr[c].re = Rand.generate_gauss();
		  psi[i].psi[mu].whr[c].im = Rand.generate_gauss();
		}
	    }
	}
    }

    void gauss(MyRand* Rand, int tid) {
      int chunk = (Z->Size)/NTHR;
      for (int i = tid*chunk; i < (tid+1)*chunk; i++) {
	for (int mu = 0; mu < dim; mu++) {
	  for (int c = 0; c < NC; c++) {
	    psi[i].psi[mu].whr[c].re = Rand[tid].generate_gauss();
	    psi[i].psi[mu].whr[c].im = Rand[tid].generate_gauss();
	  }
	}
      }
    }

    void base(int *p, int mu, int a) {
	zeros();//memset(psi, 0, sizeof(SpinColor) * Z->Size);

	psi[Z->get(p)].psi[mu].whr[a].re = 1.;
    }

    Cplx braket(int *p, int mu, int a) {
	return (psi[Z->get(p)].psi[mu].whr[a]);
    }
  
    void operator =(SpinColor_fld &scf) {
	for (int i = 0; i < Z->Size; i++)
	    for (int mu = 0; mu < dim; mu++)
		for (int a = 0; a < NC; a++)
		    psi[i].psi[mu].whr[a] = scf.psi[i].psi[mu].whr[a];
    }

    void operator +=(SpinColor_fld &scf) {
	for (int i = 0; i < Z->Size; i++)
	    for (int mu = 0; mu < dim; mu++)
		for (int a = 0; a < NC; a++)
		    psi[i].psi[mu].whr[a] += scf.psi[i].psi[mu].whr[a];
    }

    void operator -=(SpinColor_fld &scf) {
	for (int i = 0; i < Z->Size; i++)
	    psi[i] -= scf.psi[i];
    }

    void operator *=(double x) {
	for (int i = 0; i < Z->Size;i++)
	    psi[i] *= x;	
    }
    void operator *=(Cplx z) {
	for (int i = 0; i < Z->Size;i++)
	    psi[i] *= z;	
    }

    SpinColor_fld operator *(double x) {
	SpinColor_fld scf(Z);
	for (int i = 0; i < Z->Size;i++)
	    scf.psi[i] = psi[i] * x;	
	return scf;
    }

    void prout() {
	for (int i = 0; i < Z->Size; i++)
	    for (int mu = 0; mu < dim; mu++)
		for (int a = 0; a < NC; a++) {
		    if (DIFF(psi[i].psi[mu].whr[a].re, 0.0) || DIFF(psi[i].psi[mu].whr[a].im, 0.0))
			printf("%d %d %d    %+f %+f\n", i, mu, a, psi[i].psi[mu].whr[a].re, psi[i].psi[mu].whr[a].im);
		}
    }

    
    void fout(std::ofstream& file) {
      for (int i = 0; i < Z->Size; i++)
	for (int mu = 0; mu < dim; mu++)
	  for (int a = 0; a < NC; a++) {
	    if (DIFF(psi[i].psi[mu].whr[a].re, 0.0) || DIFF(psi[i].psi[mu].whr[a].im, 0.0))
	      file << i                        << " "
		   << mu                       << " "
		   << a                        << "\t"
		   << psi[i].psi[mu].whr[a].re << "\t" 
		   << psi[i].psi[mu].whr[a].im << std::endl;
	  }
    }
    
    
    SpinColor_fld *pout(dgsp &ds, eigv *ev) {
	for (int u = 0; u < ds.sz; u++)
	    psi[Z->get(ev[ds.bs + u].id.i)].psi[ev[ds.bs + u].id.mu].whr[ev[ds.bs + u].id.a] = Cplx(0., 0.);
	
	return this;
    }

    SpinColor_fld *pin(dgsp &ds, eigv *ev) {
	for (int u = 0; u < ds.bs; u++)
	    psi[Z->get(ev[u].id.i)].psi[ev[u].id.mu].whr[ev[u].id.a] = Cplx(0., 0.);

	for (int u = ds.bs + ds.sz; u < Z->Size*dim*NC; u++)
	    psi[Z->get(ev[u].id.i)].psi[ev[u].id.mu].whr[ev[u].id.a] = Cplx(0., 0.);

	return this;
    }

    SpinColor_fld *pin(dgsp &ds, eigv *ev, SpinColor_fld *star) {
	for (int u = 0; u < ds.bs; u++)
	    psi[Z->get(ev[u].id.i)].psi[ev[u].id.mu].whr[ev[u].id.a] = Cplx(0., 0.);

	for (int u = ds.bs + ds.sz; u < Z->Size*dim*NC; u++)
	    psi[Z->get(ev[u].id.i)].psi[ev[u].id.mu].whr[ev[u].id.a] = Cplx(0., 0.);

	Cplx c;
	SpinColor_fld tmp(Z);
	c = star->scal(this);
	tmp = *star;
	tmp *= c;
	*this -= tmp;
	
	return this;
    }

    void ket(int i, int mu, int a) {
	zeros();
	psi[Z->get(i)].psi[mu].whr[a].re = 1.;
    }

    void ket(indx id) {
	zeros();
	psi[Z->get(id.i)].psi[id.mu].whr[id.a].re = 1.;
    }
    
    Cplx bra(int i, int mu, int a) {
	return (psi[i].psi[mu].whr[a]);
    }

    Cplx bra(indx id) {
	return (psi[id.i].psi[id.mu].whr[id.a]);
    }

    Cplx scal(SpinColor_fld *scf) {
	Cplx s = Cplx(0, 0);
	for (int i = 0; i < Z->Size; i++)
	    for (int mu = 0; mu < dim; mu++)
		for (int a = 0; a < NC; a++)
		    s += ~(psi[i].psi[mu].whr[a]) * scf->psi[i].psi[mu].whr[a];  	
	return (s);
    }

    SpinColor near(int n, int sign, int dir) {
#ifdef PBC
      return psi[Z->near(n, sign, dir)-1];
#endif
#ifdef ABC
      int index = Z->near2(n, sign, dir);
      return SGN(index) * psi[abs(index)-1];
#endif
    }

};





