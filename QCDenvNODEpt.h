#ifndef _QCD_ENV_NODE_PT_
#define _QCD_ENV_NODE_PT_


#include <stdlib.h>
#include <stdio.h>
#include <fftw3.h>

#include "QCDenvNODE.h"

#include <Point.hpp>

#ifdef __PARALLEL_OMP__
#include <omp.h>
#endif

#define ADD 1
#define NOADD 0
#define DAG 1
#define NODAG 0


#define SGN(a) ((a) >= 0 ? +1 : -1)

const double gmT[4] = {-1.0, 1.0, -1.0, 1.0};
extern double mcpt[];



#define SWAP(a, b, s) {(s) = (a); (a) = (b); (b) = (s);}


class ptSU3_fld{
 private:
  latt *Z;
 public:
  ptSU3  *W;

  ptSU3_fld(latt* z)  {
    Z = z;
    W = new ptSU3[Z->Size];
  }
  
  ptSU3_fld(FILE*,int); // SU3_fld(&input_file,read_mode)


  ~ptSU3_fld(){
    delete[] W;
  }

  int save(char *filename) {
    FILE *filept;
    if(  (filept = fopen(filename, "wb")) == NULL)  return 1;
    fwrite(W, sizeof(ptSU3),Z->Size,filept);
    fclose(filept);
    return 0;
  }


  int load(char *filename) {
    FILE *filept;
    if( (filept = fopen(filename, "r")) == NULL)  return 1;
    fread(W, sizeof(ptSU3),Z->Size,filept);
    fclose(filept);
    return 0;
  }
  

  int load_ape(char *in){
    
    FILE *filept;
    int *xx, *YY;
    xx = new int[dim];
    YY = new int[2];
    
    char *nome;
    nome = new char[100];

    sprintf(nome,"%s.X0",in);
    if(!(filept = fopen(nome, "r"))){
      std::cout << "no configurazione " << nome << " !\n";
      exit(0);
    }
    std::cout << "File " << nome << " trovato, ora carico i dati\n";    


    
    for(YY[0] = 0; YY[0] < 8 ; YY[0]++){
      for(YY[1] = 0; YY[1] < 8 ; YY[1]++){
	
	for (int l = 0; l < 32; l++){ // t
	  for (int k = 0; k < 4; k++){ // z
	    for (int j = 0; j < 4; j++){ // y 
	      for (int i = 0; i < 8; i++){ // per noi e' x
		
		xx[3] = i;
		xx[2] = j + YY[0]*4;
		xx[1] = k + YY[1]*4;
		xx[0] = l;
		
		fread(W+Z->get(xx), sizeof(double), 2*(1+PTORD*NC*NC),filept);
		


	      }
	    }
	  }
	}
	
      }
    }

    fclose(filept);
    /***********************************/

    sprintf(nome,"%s.X1",in);
    if(!(filept = fopen(nome, "r"))){
      std::cout << "no configurazione " << nome << " !\n";
      exit(0);
    }
    std::cout << "File " << nome << " trovato, ora carico i dati\n";    


    for(YY[0] = 0; YY[0] < 8 ; YY[0]++){
      for(YY[1] = 0; YY[1] < 8 ; YY[1]++){
	
	for (int l = 0; l < 32; l++){ // t
	  for (int k = 0; k < 4; k++){ // z
	    for (int j = 0; j < 4; j++){ // y 
	      for (int i = 0; i < 8; i++){ // per noi e' x
		
		xx[3] = i+8;
		xx[2] = j + YY[0]*4;
		xx[1] = k + YY[1]*4;
		xx[0] = l;
		
		fread(W+Z->get(xx), sizeof(double), 2*(1+PTORD*NC*NC),filept);
		
	      }
	    }
	  }
	}
	
      }
    }

    fclose(filept);
    /***********************************/

    sprintf(nome,"%s.X2",in);
    if(!(filept = fopen(nome, "r"))){
      std::cout << "no configurazione " << nome << " !\n";
      exit(0);
    }
    std::cout << "File " << nome << " trovato, ora carico i dati\n";    
    

    for(YY[0] = 0; YY[0] < 8 ; YY[0]++){
      for(YY[1] = 0; YY[1] < 8 ; YY[1]++){
	
	for (int l = 0; l < 32; l++){ // t
	  for (int k = 0; k < 4; k++){ // z
	    for (int j = 0; j < 4; j++){ // y 
	      for (int i = 0; i < 8; i++){ // per noi e' x
		
		xx[3] = i+16;
		xx[2] = j + YY[0]*4;
		xx[1] = k + YY[1]*4;
		xx[0] = l;
		
		fread(W+Z->get(xx), sizeof(double), 2*(1+PTORD*NC*NC),filept);
		
	      }
	    }
	  }
	}
	
      }
    }

    fclose(filept);
    /***********************************/

    sprintf(nome,"%s.X3",in);
    if(!(filept = fopen(nome, "r"))){
      std::cout << "no configurazione " << nome << " !\n";
      exit(0);
    }
    std::cout << "File " << nome << " trovato, ora carico i dati\n";    

    for(YY[0] = 0; YY[0] < 8 ; YY[0]++){
      for(YY[1] = 0; YY[1] < 8 ; YY[1]++){
	
	for (int l = 0; l < 32; l++){ // t
	  for (int k = 0; k < 4; k++){ // z
	    for (int j = 0; j < 4; j++){ // y 
	      for (int i = 0; i < 8; i++){ // per noi e' x
		
		xx[3] = i+24;
		xx[2] = j + YY[0]*4;
		xx[1] = k + YY[1]*4;
		xx[0] = l;
		
		fread(W+Z->get(xx), sizeof(double), 2*(1+PTORD*NC*NC),filept);
		
	      }
	    }
	  }
	}
	
      }
    }

    fclose(filept);
    delete [] nome;

    return 1;
  }




  ptSU3* handle(){ return W; }

  ptSU3 get(int *);

  ptSU3 get(int n){
    return W[Z->L[n][4]];
  }
  
  ptSU3 get(int n,int step,int dir){
    return W[Z->get(n, step, dir)];
  }


  friend int get(ptSU3_fld *, int*);
  friend int get(ptSU3_fld *, int, int, int);
  friend int get(ptSU3_fld *, int);

};
  
inline ptSU3 ptSU3_fld::get(int *n) {
    return W[Z->get(n)];
}

inline int get(ptSU3_fld *W, int *n) {
    return (W->Z)->get(n);

}

inline int get(ptSU3_fld *W,int n,int step,int dir){
  return (W->Z)->get(n,step,dir);
}

inline int get(ptSU3_fld *W,int n){
  return (W->Z)->get(n);
}

// ------------- end class pt SU3_fld  --------------



class ptGluon_fld{
 private:


 public:
  latt   *Z;
  ptGluon  *W;

  ptGluon_fld(latt *z){
    Z = z;
    W = new ptGluon[z->Size];
  }
  
  ~ptGluon_fld(){
    delete [] W;
  }

  ptGluon_fld(FILE*,int); // Gluon_fld(&input_file,read_mode)
  
  int save(char *filename) {
    FILE *filept;	
    if( (filept = fopen(filename, "w")) == NULL)  return 1;

    fwrite(W,sizeof(ptGluon), Z->Size, filept);
    fclose(filept); 

    return 0;
  }

  int load(char *filename) {
    FILE *filept;
    if( (filept = fopen(filename, "r")) == NULL)  {
      cout << filept << endl;
      return 1;
    }
    fread(W,sizeof(ptGluon), Z->Size,filept);

    fclose(filept);
    return 0;
  }

  // salva anche la placchetta, in modo
  // da avere una verifica di correttezza
  // sulla configurazione
  int save(char *filename, Cplx *placchetta) {
    FILE *filept;	
    if( (filept = fopen(filename, "w")) == NULL)  return 1;
    fwrite(W,sizeof(ptGluon), Z->Size, filept);
    fwrite(placchetta, sizeof(Cplx), PTORD+1 , filept);
    
    fclose(filept);
    return 0;
  }

  // legge anche la placchetta, in modo
  // da avere una verifica di correttezza
  // sulla configurazione
  int load(char *filename, Cplx *placchetta) {
    FILE *filept;
    if( (filept = fopen(filename, "r") ) == NULL )
      {
	return 1;
      }

    fread(W,sizeof(ptGluon), Z->Size,filept); 
    fread(placchetta, sizeof(Cplx), PTORD+1 , filept);

    fclose(filept);
    return 0;
  }


  int load(char *filename, Cplx *placchetta, int rank) {
    FILE *filept;
    if( (filept = fopen(filename, "r") ) == NULL )
      {
	return 1;
      }

    long offset = Z->Sz[0]*rank*sizeof(ptGluon);
    fseek(filept, offset, SEEK_SET);
    fread(W,sizeof(ptGluon), Z->Size,filept); 

    offset = -(PTORD+1)*sizeof(Cplx);
    fseek(filept, offset, SEEK_END);
    
    fread(placchetta, sizeof(Cplx), PTORD+1 , filept);

    fclose(filept);
    return 0;
  }
  


  int load_ape(char* in){
    
    FILE *filept;
    int *xx, *YY;
    xx = new int[dim];
    YY = new int[2];
    
    char *nome;
    nome = new char[100];
    sprintf(nome,"%s.X0",in);

    if(!(filept = fopen(nome, "r"))){
      std::cout << "no configurazione " << nome << " !\n";
      exit(0);
    }
    std::cout << "File " << nome << " trovato, ora carico i dati\n";
    
    for(YY[0] = 0; YY[0] < 8 ; YY[0]++){
      for(YY[1] = 0; YY[1] < 8 ; YY[1]++){
	
	for (int l = 0; l < 32; l++){ // t
	  for (int k = 0; k < 4; k++){ // z
	    for (int j = 0; j < 4; j++){ // y 
	      for (int i = 0; i < 8; i++){ // per noi e' x
		
		xx[3] = l;
		xx[2] = k + YY[1]*4;
		xx[1] = j + YY[0]*4;
		xx[0] = i;
		
		fread(W+Z->get(xx), sizeof(double), 2*(1+PTORD*NC*NC)*4,filept);
		


	      }
	    }
	  }
	}
	
      }
    }

    fclose(filept);
    /***********************************/

    sprintf(nome,"%s.X1",in);
    if(!(filept = fopen(nome, "r"))){
      std::cout << "no configurazione " << nome << " !\n";
      exit(0);
    }
    std::cout << "File " << nome << " trovato, ora carico i dati\n";
    

    for(YY[0] = 0; YY[0] < 8 ; YY[0]++){
      for(YY[1] = 0; YY[1] < 8 ; YY[1]++){
	
	for (int l = 0; l < 32; l++){ // t
	  for (int k = 0; k < 4; k++){ // z
	    for (int j = 0; j < 4; j++){ // y 
	      for (int i = 0; i < 8; i++){ // per noi e' x
		
		xx[3] = l;
		xx[2] = k + YY[1]*4;
		xx[1] = j + YY[0]*4;
		xx[0] = i+8;
/* 		xx[3] = l; */
/* 		xx[2] = k + YY[0]*4; */
/* 		xx[1] = j + YY[1]*4; */
/* 		xx[0] = i+8; */
		
		fread(W+Z->get(xx), sizeof(double), 2*(1+PTORD*NC*NC)*4,filept);
		
	      }
	    }
	  }
	}
	
      }
    }

    fclose(filept);
    /***********************************/

    sprintf(nome,"%s.X2",in);
    if(!(filept = fopen(nome, "r"))){
      std::cout << "no configurazione " << nome << " !\n";
      exit(0);
    }
    std::cout << "File " << nome << " trovato, ora carico i dati\n";
    

    for(YY[0] = 0; YY[0] < 8 ; YY[0]++){
      for(YY[1] = 0; YY[1] < 8 ; YY[1]++){
	
	for (int l = 0; l < 32; l++){ // t
	  for (int k = 0; k < 4; k++){ // z
	    for (int j = 0; j < 4; j++){ // y 
	      for (int i = 0; i < 8; i++){ // per noi e' x
		
		xx[3] = l;
		xx[2] = k + YY[1]*4;
		xx[1] = j + YY[0]*4;
		xx[0] = i+16;
/* 		xx[3] = l; */
/* 		xx[2] = k + YY[0]*4; */
/* 		xx[1] = j + YY[1]*4; */
/* 		xx[0] = i+16; */
		
		fread(W+Z->get(xx), sizeof(double), 2*(1+PTORD*NC*NC)*4,filept);
		
	      }
	    }
	  }
	}
	
      }
    }

    fclose(filept);
    /***********************************/

    sprintf(nome,"%s.X3",in);
    if(!(filept = fopen(nome, "r"))){
      std::cout << "no configurazione " << nome << " !\n";
      exit(0);
    }
    std::cout << "File " << nome << " trovato, ora carico i dati\n";    

    for(YY[0] = 0; YY[0] < 8 ; YY[0]++){
      for(YY[1] = 0; YY[1] < 8 ; YY[1]++){
	
	for (int l = 0; l < 32; l++){ // t
	  for (int k = 0; k < 4; k++){ // z
	    for (int j = 0; j < 4; j++){ // y 
	      for (int i = 0; i < 8; i++){ // per noi e' x
		
		xx[3] = l;
		xx[2] = k + YY[1]*4;
		xx[1] = j + YY[0]*4;
		xx[0] = i+24;
/* 		xx[3] = l; */
/* 		xx[2] = k + YY[0]*4; */
/* 		xx[1] = j + YY[1]*4; */
/* 		xx[0] = i+24; */
		
		fread(W+Z->get(xx), sizeof(double), 2*(1+PTORD*NC*NC)*4,filept);
		
	      }
	    }
	  }
	}
	
      }
    }

    fclose(filept);

    delete [] nome;
    /*
    ptSU3 scambio;

    for(int sito = 0; sito < Z->Size; sito++){
      scambio      = W[sito].U[0];
      W[sito].U[0] = W[sito].U[3];
      W[sito].U[3] = scambio;
      scambio      = W[sito].U[1];
      W[sito].U[1] = W[sito].U[2];
      W[sito].U[2] = scambio;
    }
    */
    return 1;
  }


  void prout(){
    for(int i = 0; i < Z->Size; i++){
      printf("i = %d\n",i);
      
      //W[i].prout();

      printf("\n\n");
    }
  }


  void operator=(const ptGluon_fld& A) {
    Z = A.Z;
    for(int i = 0; i < Z->Size; i++)
      W[i] = A.W[i];
  }

  ptSU3* handle(){ return (ptSU3*)(&(W[0])); }
  
  ptGluon get(int*);

  ptGluon get(int n, int step, int dir) {return W[Z->get(n, step, dir)];}

  ptSU3 get(int*, int);

  ptSU3 get(int n, int mu){
      return W[Z->get(n)][mu];
  }

  ptSU3 get(int n, int step, int dir, int mu){
    return W[Z->get(n, step, dir)][mu];
  }
  
  friend int get(ptGluon_fld *W, int n, int mu) {
    return ( 4*((W->Z)->get(n))+mu );
  }

  friend int get(ptGluon_fld *, int *);

  friend int get(ptGluon_fld *W, int n) {
    return ( 4*((W->Z)->get(n)) );
  }

  friend int get(ptGluon_fld *W, int n, int step , int dir, int mu) {
    return (4*((W->Z)->get(n,step,dir))+mu);
  }

  friend int get(ptGluon_fld *W, int n, int step1 , int dir1,
		 int step2, int dir2, int mu) {
    return (4*((W->Z)->get( ((W->Z)->get(n, step1, dir1)), step2, dir2))+mu);
  }


  ptGluon near(int n, int sign, int dir) {
    return W[Z->near(n, sign, dir)-1];
  }    

  ptSU3 staple(int*, int, int); //(x[dim], mu, nu)
  
  inline ptSU3 staple(int, int, int);  //(n, mu, nu)
  inline ptSU3 staple2x1(int, int, int); //(n, mu, nu)
  //  inline void wilson2x2(ptSU3&, int, int, int); //(n, mu, nu)

  ptSU3& operator()(const pt::Point&, const pt::Direction&); 
  const ptSU3& operator()(const pt::Point&, const pt::Direction&) const;
  pt::Point mk_point(int t, int x, int y, int z) const;

}; // class ptGluon_fld


// need to implement out of class
// otherwise, we would get an "incomplete type" error

inline ptSU3& ptGluon_fld::operator()(const pt::Point& n, const pt::Direction& mu) {
  //return mu.deref_bkw(n.deref(W));
  return mu.deref_bkw<ptSU3&, ptGluon&>(n.deref<ptGluon>(W));
}

inline pt::Point ptGluon_fld::mk_point(int t, int x, int y, int z) const {
  int a[4] = {t,x,y,z};
  return pt::Point(Z->get(a), Z->L.begin());
}

inline const ptSU3& ptGluon_fld::operator()
  (const pt::Point& n, const pt::Direction& mu) const {
  return mu.deref_bkw<const ptSU3&, const ptGluon&>(n.deref(W));
}

inline ptGluon ptGluon_fld::get(int *n) {
    return W[Z->get(n)];
  }
  
inline ptSU3 ptGluon_fld::get(int *n, int mu) {
    return W[Z->get(n)][mu];
  }
  
inline int get(ptGluon_fld *W, int *n) {
    return W->Z->get(n);
}

ptSU3 ptGluon_fld::staple(int n, int mu, int nu){

  /*

    b ->
    ^    ^
    |    |
    n == a

    n == 
    ^    ^
    |    |
    b -> a

   */


  return( W[Z->get(n, 1, mu)][nu]*
	  dag(W[Z->L[n][4]][nu]*
	      W[Z->get(n, 1, nu)][mu]) 
	  +
	  dag(W[Z->get(n, -1, nu)][mu]*
	      W[Z->get(n, 1, mu, -1, nu)][nu])*
	  W[Z->get(n, -1, nu)][nu]);
}

ptSU3 ptGluon_fld::staple2x1(int n, int mu, int nu){

  int curr = Z->L[n][4];
  int a,b,c,d,e,f,g,h;
  ptSU3 tmp;

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
  g = Z->L[f][5+nu];

  tmp = ( W[a][nu] * W[b][nu] * dag(W[curr][nu] * W[c][nu] * W[d][mu] ) +
          W[a][mu] * W[e][nu] * dag(W[curr][nu] * W[c][mu] * W[b][mu] ) +
          W[a][nu] * dag(W[f][nu] * W[g][mu] * W[c][mu]) * W[f][mu] );
  
  //  a = Z->L[curr][5+mu];
  b = Z->L[a][nu  ];
  c = Z->L[curr][nu  ];
  d = Z->L[c][nu  ];
  e = Z->L[b][nu  ];
  f = Z->L[b][5+mu];
  g = Z->L[c][mu  ];
  h = Z->L[curr][mu   ];

  /*   b^[nu] e^[nu] d^[mu] d[nu] c[nu] = (d[mu] e[nu] b[nu] )^ d[nu] c[nu] */
  /*   a[mu] f^[nu] b^[mu] c^[mu] c[nu] = a[mu] (c[mu] b[mu] f[nu] )^ c[nu] */
  /*   b^[nu] c^[mu] g^[mu] g[nu] h[mu] = (g[mu] c[mu] b[nu] )^ g[nu] h[mu] */
    
  tmp += ( dag(W[d][mu] * W[e][nu] * W[b][nu]) * W[d][nu] * W[c][nu]  +
           W[a][mu] * dag(W[c][mu] * W[b][mu] * W[f][nu]) * W[c][nu]  +
           dag(W[g][mu] * W[c][mu] * W[b][nu]) * W[g][nu] * W[h][mu] );

  return tmp;
}


/* void ptGluon_fld::wilson2x2(ptSU3& tmp, int n, int mu, int nu){ */
  
/*   //  int curr = Z->L[n][4]; */
/*   int a,b,c,d,e,f; */
/*   //  ptSU3 tmp; */
  
/*   /\* Wilson Loop 2x2. */
     
/*      e -> f ->  */
/*      ^         ^ */
/*      |         | */
/*      d         c */
/*      ^         ^ */
/*      |         | */
/*      n -> a -> b */
     
/*   *\/ */
  
/*   a = Z->L[n][5+mu]; */
/*   b = Z->L[a][5+mu]; */
/*   c = Z->L[b][5+nu]; */
/*   d = Z->L[n][5+nu]; */
/*   e = Z->L[d][5+nu]; */
/*   f = Z->L[e][5+mu]; */


/*   tmp += ( W[n][mu] * W[a][mu] * W[b][nu] * W[c][nu] * */
/* 	   dag( W[n][nu] * W[d][nu] * W[e][mu] * W[f][mu]) ); */
  
  
/*   //  return tmp; */
/* } */



// ********* end class Gluon_fld **************

class ptSpinColor_fld {
  
public:
  latt *Z;
  ptSpinColor *psi;

  const ptSpinColor& operator[](int i) const { return psi[i]; };
  ptSpinColor& operator[](int i) { return psi[i]; };
  
  fftw_plan plan[2];
  SpinColor_fld *scfld;

  ptSpinColor_fld(latt *z) {
    Z = z;
    
    psi = new ptSpinColor[Z->Size];
    scfld = new SpinColor_fld(Z);
    
#ifdef __PARALLEL_OMP__
    if( !fftw_init_threads() ) 
      {
	printf("\n");
	printf("######################################################################\n");
	printf("There was an error in multi-threaded FFTW, program will be terminated.\n");
	printf("######################################################################\n");
	printf("\n");

      }
    else{
      fftw_plan_with_nthreads(NTHR); 
    }
#endif

    plan[0] = fftw_plan_many_dft(dim, Z->Sz, NC*dim,
				 (fftw_complex *) scfld->psi, NULL, NC*dim, 1,
				 (fftw_complex *) scfld->psi, NULL, NC*dim, 1,
				 FFTW_FORWARD, FFTW_MEASURE);
    plan[1] = fftw_plan_many_dft(dim, Z->Sz, NC*dim,
				 (fftw_complex *) scfld->psi, NULL, NC*dim, 1,
				 (fftw_complex *) scfld->psi, NULL, NC*dim, 1,
				 FFTW_BACKWARD, FFTW_MEASURE);
  }
  
  
  ~ptSpinColor_fld() {
    delete [] psi;
    delete scfld;
#ifndef __PARALLEL_OMP__
    fftw_cleanup();
#else
    fftw_cleanup_threads();
#endif
  }
  

  
  int save(char *filename) {
    FILE *filept;	
    if( (filept = fopen(filename, "w")) == NULL) {
      return 1;
    }
    if( fwrite(psi,sizeof(ptSpinColor), Z->Size,filept) ) {
      fclose(filept);
      return 1;
    }
    fclose(filept);
    return 0;
  }
  
  int load(char *filename) {
    FILE *filept;
    if( (filept = fopen(filename, "r")) == NULL) return 1;
    if( fread(psi,sizeof(ptSpinColor), Z->Size,filept) ) {
      fclose(filept);
      return 1;
    }
    fclose(filept);
    return 0;
  }
  
  
  void prout(){
    for(int i = 0; i < Z->Size; i++){
      printf("i = %d\n",i);
      for(int mu = 0; mu < dim; mu++){
	printf("\tmu = %d\n",mu);
	for(int ord = 0; ord < PTORD+1; ord++){
	  printf("\t\tord = %d\n",ord);
	  psi[i][mu][ord].prout();
	  printf("\n");
	}
	printf("\n");
      }
      printf("\n\n");
    }
  }
  
  void prout(int o){
    for(int i = 0; i < Z->Size; i++){
      for(int mu = 0; mu < dim; mu++){
	for( int a = 0; a < NC; a++){
	  if(fabs(psi[i][mu][o].whr[a].re)+
	     fabs(psi[i][mu][o].whr[a].im) > 1e-10 ){
	    printf("%d\t%d\t%d:\t",i,mu,a);
	    psi[i][mu][o].whr[a].prout();
	    printf("\n");
	  }
	}
      }
    }
  }
    


  void fout(std::ofstream& file) {
    for(int ord = 0; ord < PTORD+1; ord++){
      file << "ord = " << ord << endl << endl;
      for (int i = 0; i < Z->Size; i++)
	for (int mu = 0; mu < dim; mu++)
	  for (int a = 0; a < NC; a++) {
	    if(fabs(psi[i][mu][ord].whr[a].re)+
	       fabs(psi[i][mu][ord].whr[a].im) > 1e-10 )
	      file << i                        << " "
		   << mu                       << " "
		   << a                        << "\t"
		   << psi[i][mu][ord].whr[a].re << "\t" 
		   << psi[i][mu][ord].whr[a].im << std::endl;
	  }
      file << endl;
    }
  }
  

  void fout(std::ofstream& file, int ord) {
      for (int i = 0; i < Z->Size; i++)
	for (int mu = 0; mu < dim; mu++)
	  for (int a = 0; a < NC; a++) {
	    if(fabs(psi[i][mu][ord].whr[a].re)+
	       fabs(psi[i][mu][ord].whr[a].im) > 1e-10 )
	      file << i                        << " "
		   << mu                       << " "
		   << a                        << "\t"
		   << psi[i][mu][ord].whr[a].re << "\t" 
		   << psi[i][mu][ord].whr[a].im << std::endl;
	  }
  }





  void operator=(ptSpinColor_fld &ptscfld) {
    Z = ptscfld.Z;
    for(int i = 0; i < Z->Size; i++)
      psi[i] = ptscfld[i];
  }
  

  void scalToPt(int ord) {
    
    for(int i = 0; i < Z->Size; i++){
      for(int mu = 0; mu < dim; mu++){
        psi[i][mu][ord] = (*scfld)[i][mu];
      }
    }

  }


  ptCVector* handle() {
    return (ptCVector *) (& (psi[0]));
  }
  
  ptSpinColor get(int *n){ return psi[Z->get(n)]; }
  
  ptSpinColor get(int n, int step, int dir) {
    return psi[Z->get(n, step, dir)];
  }
  
  ptSpinColor get(int n) {
    return psi[Z->get(n)];
  }
  
  ptCVector get(int n, int dir) {
    return psi[Z->get(n)][dir];
  }


  ptSpinColor near(int n, int sign, int dir) {
#ifdef PBC
    return psi[Z->near(n, sign, dir)-1];
#endif
#ifdef ABC
    int index = Z->near2(n, sign, dir);
    return SGN(index) * psi[abs(index)-1];
#endif
  }
  


  void fft(int index) {
    
#ifdef ABC
    int coord[4];
    switch(index){
      
    case 0:
      for (int i = 0; i < Z->Size; i++){
	Z->get(i,coord);
	scfld->psi[i] *= Cplx(cos(M_PI*coord[0]/Z->Sz[0]), sin(M_PI*coord[0]/Z->Sz[0]));
      }
      fftw_execute(plan[0]);
      break;
      
    case 1:
      fftw_execute(plan[1]);
      for (int i = 0; i < Z->Size; i++){
	for (int mu = 0; mu < dim; mu++){
	  for (int a = 0; a < NC; a++) {
	    scfld->psi[i][mu].whr[a].re /= (double) Z->Size;
	    scfld->psi[i][mu].whr[a].im /= (double) Z->Size;
	  }
	}
	
	Z->get(i,coord);
	scfld->psi[i] *= Cplx(cos(M_PI*coord[0]/Z->Sz[0]), -sin(M_PI*coord[0]/Z->Sz[0]));
      }
      break;
    }

#elif defined PBC
    fftw_execute(plan[index]);
    if(index==1){
      for (int i = 0; i < Z->Size; i++){
	for (int mu = 0; mu < dim; mu++){
	  for (int a = 0; a < NC; a++) {
	    (*scfld)[i][mu][a] /= (double) Z->Size;
	  }
	}
      }
    }
#endif

  }




  void fftT(int index) {
    
#ifdef ABC
    int coord[4];
    switch(index){
      
    case 0:
      for (int i = 0; i < Z->Size; i++){
	Z->get(i,coord);
	scfld->psi[i] *= Cplx(cos(M_PI*coord[0]/Z->Sz[0]), sin(M_PI*coord[0]/Z->Sz[0]));
      }
      fftw_execute(plan[1]);
      break;
      
    case 1:
      fftw_execute(plan[0]);
      for (int i = 0; i < Z->Size; i++){
	for (int mu = 0; mu < dim; mu++){
	  for (int a = 0; a < NC; a++) {
	    scfld->psi[i][mu].whr[a].re /= (double) Z->Size;
	    scfld->psi[i][mu].whr[a].im /= (double) Z->Size;
	  }
	}
	
	Z->get(i,coord);
	scfld->psi[i] *= Cplx(cos(M_PI*coord[0]/Z->Sz[0]), -sin(M_PI*coord[0]/Z->Sz[0]));
      }
      break;
    }
#elif defined PBC
    fftw_execute(plan[(index+1)%2]);
    if(index==1){
      for (int i = 0; i < Z->Size; i++){
	for (int mu = 0; mu < dim; mu++){
	  for (int a = 0; a < NC; a++) {
	    scfld->psi[i][mu].whr[a] /= (double) Z->Size;
	  }
	}
      }
    }
#endif
  }
  

#ifndef __PARALLEL_OMP__
  
  void M0inv(){
    double diag, den;
    Cplx im_(0, -1);
    int i = 0, pp[4];

#ifdef PBC
    scfld->psi[Z->get(0)] *= Cplx(0.0,0.0);
#endif
    for (pp[0] = 0; pp[0] < Z->Sz[0]; pp[0]++) {
      for (pp[1] = 0; pp[1] < Z->Sz[1]; pp[1]++) {
	for (pp[2] = 0; pp[2] < Z->Sz[2]; pp[2]++) {		
	  for (pp[3] = 0; pp[3] < Z->Sz[3]; pp[3]++) {
#ifdef PBC
	    if( (pp[0]+pp[1]+pp[2]+pp[3])!= 0){
#endif
	      i = Z->get(pp);
	      diag = MBARE + 0.5*(Z->p2hat[0][pp[0]] + Z->p2hat[1][pp[1]] + 
				  Z->p2hat[2][pp[2]] + Z->p2hat[3][pp[3]]);
	      den = diag*diag + (Z->p2bar[0][pp[0]] + Z->p2bar[1][pp[1]] + 
				 Z->p2bar[2][pp[2]] + Z->p2bar[3][pp[3]]);
	      
	      scfld->psi[i] = ( ((scfld->psi[i]).gmleft(0) * (Z->pbar[0][pp[0]]) +
				 (scfld->psi[i]).gmleft(1) * (Z->pbar[1][pp[1]]) +
				 (scfld->psi[i]).gmleft(2) * (Z->pbar[2][pp[2]]) +
				 (scfld->psi[i]).gmleft(3) * (Z->pbar[3][pp[3]]))
				*im_+ (scfld->psi[i])*diag) / den;
#ifdef PBC
	    }
#endif
	  }
	}
      }
    }
  }
  
#else

  void M0inv(){
    double diag, den;
    Cplx im_(0, -1);
    int i = 0, pp[NTHR][4],tid;
    int chunk = (Z->Size)/NTHR;
#ifdef ABC

#pragma omp parallel num_threads(NTHR) private(tid,i, diag, den)
    {
      tid = omp_get_thread_num();
      for(i = tid*chunk; i < (tid+1)*chunk; i++){
	Z->get(i,pp[tid]);
	diag = MBARE + 0.5*(Z->p2hat[0][pp[tid][0]] + Z->p2hat[1][pp[tid][1]] + 
			    Z->p2hat[2][pp[tid][2]] + Z->p2hat[3][pp[tid][3]]);
	den = diag*diag + (Z->p2bar[0][pp[tid][0]] + Z->p2bar[1][pp[tid][1]] + 
			 Z->p2bar[2][pp[tid][2]] + Z->p2bar[3][pp[tid][3]]);
	
	scfld->psi[i] = ( ((scfld->psi[i]).gmleft(0) * (Z->pbar[0][pp[tid][0]]) +
			   (scfld->psi[i]).gmleft(1) * (Z->pbar[1][pp[tid][1]]) +
			   (scfld->psi[i]).gmleft(2) * (Z->pbar[2][pp[tid][2]]) +
			   (scfld->psi[i]).gmleft(3) * (Z->pbar[3][pp[tid][3]]))
			  *im_+ (scfld->psi[i])*diag) / den;
      } // siti
    } //parallel
#elif defined PBC
    scfld->psi[Z->get(0)] *= Cplx(0.0,0.0);

#pragma omp parallel num_threads(NTHR) private(tid,i, diag, den)
    {
      tid = omp_get_thread_num();
      for(i = tid*chunk; i < (tid+1)*chunk; i++){
	Z->get(i,pp[tid]);
	if( (pp[tid][0]+pp[tid][1]+pp[tid][2]+pp[tid][3])!= 0){
	  i = Z->get(pp[tid]);
	  diag = MBARE + 0.5*(Z->p2hat[0][pp[tid][0]] + Z->p2hat[1][pp[tid][1]] + 
			      Z->p2hat[2][pp[tid][2]] + Z->p2hat[3][pp[tid][3]]);
	  den = diag*diag + (Z->p2bar[0][pp[tid][0]] + Z->p2bar[1][pp[tid][1]] + 
			     Z->p2bar[2][pp[tid][2]] + Z->p2bar[3][pp[tid][3]]);
	  
	  scfld->psi[i] = ( ((scfld->psi[i]).gmleft(0) * (Z->pbar[0][pp[tid][0]]) +
			     (scfld->psi[i]).gmleft(1) * (Z->pbar[1][pp[tid][1]]) +
			     (scfld->psi[i]).gmleft(2) * (Z->pbar[2][pp[tid][2]]) +
			     (scfld->psi[i]).gmleft(3) * (Z->pbar[3][pp[tid][3]]))
			    *im_+ (scfld->psi[i])*diag) / den;
	}
      } // fine siti
    } // end parallel
#endif
  }
    
#endif
  
  
  
  void M0invT(){
  
    double diag, den;
    Cplx im_(0, -1);
    int i = 0, pp[4];
  
#ifdef ABC    
    for (pp[0] = 0; pp[0] < Z->Sz[0]; pp[0]++) {
      for (pp[1] = 0; pp[1] < Z->Sz[1]; pp[1]++) {
	for (pp[2] = 0; pp[2] < Z->Sz[2]; pp[2]++) {		
	  for (pp[3] = 0; pp[3] < Z->Sz[3]; pp[3]++) {
	    i = Z->get(pp);
	    diag = MBARE + 0.5*(Z->p2hat[0][pp[0]] + Z->p2hat[1][pp[1]] + 
				Z->p2hat[2][pp[2]] + Z->p2hat[3][pp[3]]);
	    den = diag*diag + (Z->p2bar[0][pp[0]] + Z->p2bar[1][pp[1]] + 
			       Z->p2bar[2][pp[2]] + Z->p2bar[3][pp[3]]);
	    
	    scfld->psi[i] = ( ((scfld->psi[i]).gmleft(0)*gmT[0] * (Z->pbar[0][pp[0]]) +
			       (scfld->psi[i]).gmleft(1)*gmT[1] * (Z->pbar[1][pp[1]]) +
			       (scfld->psi[i]).gmleft(2)*gmT[2] * (Z->pbar[2][pp[2]]) +
			       (scfld->psi[i]).gmleft(3)*gmT[3] * (Z->pbar[3][pp[3]]))
			      *im_+ (scfld->psi[i])*diag) / den;
	  }
	}
      }
    }
#elif defined PBC
    scfld->psi[Z->get(0)] *= Cplx(0.0,0.0);
    for (pp[0] = 0; pp[0] < Z->Sz[0]; pp[0]++) {
      for (pp[1] = 0; pp[1] < Z->Sz[1]; pp[1]++) {
	for (pp[2] = 0; pp[2] < Z->Sz[2]; pp[2]++) {		
	  for (pp[3] = 0; pp[3] < Z->Sz[3]; pp[3]++) {
	    if( (pp[0]+pp[1]+pp[2]+pp[3])!= 0){
	      i = Z->get(pp);
	      diag = MBARE + 0.5*(Z->p2hat[0][pp[0]] + Z->p2hat[1][pp[1]] + 
				  Z->p2hat[2][pp[2]] + Z->p2hat[3][pp[3]]);
	      den = diag*diag + (Z->p2bar[0][pp[0]] + Z->p2bar[1][pp[1]] + 
				 Z->p2bar[2][pp[2]] + Z->p2bar[3][pp[3]]);
	      
	      scfld->psi[i] = ( ((scfld->psi[i]).gmleft(0)*gmT[0] * (Z->pbar[0][pp[0]]) +
				 (scfld->psi[i]).gmleft(1)*gmT[1] * (Z->pbar[1][pp[1]]) +
				 (scfld->psi[i]).gmleft(2)*gmT[2] * (Z->pbar[2][pp[2]]) +
				 (scfld->psi[i]).gmleft(3)*gmT[3] * (Z->pbar[3][pp[3]]))
				*im_+ (scfld->psi[i])*diag) / den;
	    }
	  }
	}
      }
    }
#endif
  }








#ifndef __PARALLEL_OMP__

 void fillPT(int ptmax, ptGluon_fld &Umu){
   
   fft(0);
   M0inv();
   fft(1);
   
   int site_c, point_up, point_dn;
   SpinColor Xi1, Xi2;

   // per copiare non mi interessa l'ordine in cui accedo ai siti
   for(int i = 0; i < Z->Size; i++){
     for(int mu = 0; mu < dim; mu++){
       psi[i][mu][0] = scfld->psi[i][mu];
     }
   }

   for(int jord = 1; jord <= ptmax; jord++){
      
     for(int curr = 0; curr < Z->Size; curr++){
	
       site_c = Z->get(curr);
	      
       // moltiplicazione perturbativa per M
       for( int kord = 0; kord < jord; kord++) {
	  
	 // massa critica
	 for (int mu = 0; mu < dim; mu++ ) {		  
	   psi[site_c][mu][jord] += mcpt[jord-kord]*psi[site_c][mu][kord];
	 }

#ifndef _USE_HALFSPINOR_
	 for( int mu = 0; mu < dim; mu++){
	    
	   point_up = Umu.Z->get(curr, 1, mu);
	   point_dn = Umu.Z->get(curr,-1, mu);
	    
	   psi[point_dn].uno_p_gmu(Xi1,mu,kord);
	   psi[point_up].uno_m_gmu(Xi2,mu,kord);
	    
	   for( int nu = 0; nu < dim; nu++ ){
	      
	     psi[site_c][nu][jord] -= .5*dag(Umu.W[point_dn][mu][jord-kord-1])*Xi1[nu];
	     psi[site_c][nu][jord] -= .5*( Umu.W[site_c][mu][jord-kord-1] )*Xi2[nu];
	   } //nu
	    
	 } // mu
#else	    
	 // -- Half Spinor --
	 // see arXiv:0905.3331v1 [hep-lat]

	 // mu = 0
	 point_up = Umu.Z->get(curr, 1, 0);
	 point_dn = Umu.Z->get(curr,-1, 0);
	 
	 Xi1[0].whr[0].re = psi[point_dn][0][kord].whr[0].re + psi[point_dn][3][kord].whr[0].im;    Xi1[0].whr[0].im = psi[point_dn][0][kord].whr[0].im - psi[point_dn][3][kord].whr[0].re;
	 Xi1[0].whr[1].re = psi[point_dn][0][kord].whr[1].re + psi[point_dn][3][kord].whr[1].im;    Xi1[0].whr[1].im = psi[point_dn][0][kord].whr[1].im - psi[point_dn][3][kord].whr[1].re;
	 Xi1[0].whr[2].re = psi[point_dn][0][kord].whr[2].re + psi[point_dn][3][kord].whr[2].im;    Xi1[0].whr[2].im = psi[point_dn][0][kord].whr[2].im - psi[point_dn][3][kord].whr[2].re;
	 Xi1[1].whr[0].re = psi[point_dn][1][kord].whr[0].re + psi[point_dn][2][kord].whr[0].im;    Xi1[1].whr[0].im = psi[point_dn][1][kord].whr[0].im - psi[point_dn][2][kord].whr[0].re;
	 Xi1[1].whr[1].re = psi[point_dn][1][kord].whr[1].re + psi[point_dn][2][kord].whr[1].im;    Xi1[1].whr[1].im = psi[point_dn][1][kord].whr[1].im - psi[point_dn][2][kord].whr[1].re;
	 Xi1[1].whr[2].re = psi[point_dn][1][kord].whr[2].re + psi[point_dn][2][kord].whr[2].im;    Xi1[1].whr[2].im = psi[point_dn][1][kord].whr[2].im - psi[point_dn][2][kord].whr[2].re;

	 Xi1[0] = dag(Umu.W[point_dn][0][jord-kord-1])*Xi1[0];
	 Xi1[1] = dag(Umu.W[point_dn][0][jord-kord-1])*Xi1[1];

	 Xi2[0].whr[0].re = psi[point_up][0][kord].whr[0].re - psi[point_up][3][kord].whr[0].im;      Xi2[0].whr[0].im = psi[point_up][0][kord].whr[0].im + psi[point_up][3][kord].whr[0].re;
	 Xi2[0].whr[1].re = psi[point_up][0][kord].whr[1].re - psi[point_up][3][kord].whr[1].im;      Xi2[0].whr[1].im = psi[point_up][0][kord].whr[1].im + psi[point_up][3][kord].whr[1].re;
	 Xi2[0].whr[2].re = psi[point_up][0][kord].whr[2].re - psi[point_up][3][kord].whr[2].im;      Xi2[0].whr[2].im = psi[point_up][0][kord].whr[2].im + psi[point_up][3][kord].whr[2].re;
	 Xi2[1].whr[0].re = psi[point_up][1][kord].whr[0].re - psi[point_up][2][kord].whr[0].im;      Xi2[1].whr[0].im = psi[point_up][1][kord].whr[0].im + psi[point_up][2][kord].whr[0].re;
	 Xi2[1].whr[1].re = psi[point_up][1][kord].whr[1].re - psi[point_up][2][kord].whr[1].im;      Xi2[1].whr[1].im = psi[point_up][1][kord].whr[1].im + psi[point_up][2][kord].whr[1].re;
	 Xi2[1].whr[2].re = psi[point_up][1][kord].whr[2].re - psi[point_up][2][kord].whr[2].im;      Xi2[1].whr[2].im = psi[point_up][1][kord].whr[2].im + psi[point_up][2][kord].whr[2].re;

	 Xi2[0] = Umu.W[curr][0][jord-kord-1]*Xi2[0];
	 Xi2[1] = Umu.W[curr][0][jord-kord-1]*Xi2[1];
	 
	 // reconstruct
	 psi[site_c][0][jord] -= .5*(Xi1[0] + Xi2[0]);
	 psi[site_c][1][jord] -= .5*(Xi1[1] + Xi2[1]);

	 psi[site_c][2][jord].whr[0].re += .5*( Xi1[1].whr[0].im - Xi2[1].whr[0].im );
	 psi[site_c][2][jord].whr[0].im -= .5*( Xi1[1].whr[0].re - Xi2[1].whr[0].re );
	 psi[site_c][2][jord].whr[1].re += .5*( Xi1[1].whr[1].im - Xi2[1].whr[1].im );
	 psi[site_c][2][jord].whr[1].im -= .5*( Xi1[1].whr[1].re - Xi2[1].whr[1].re );
	 psi[site_c][2][jord].whr[2].re += .5*( Xi1[1].whr[2].im - Xi2[1].whr[2].im );
	 psi[site_c][2][jord].whr[2].im -= .5*( Xi1[1].whr[2].re - Xi2[1].whr[2].re );

	 psi[site_c][3][jord].whr[0].re += .5*( Xi1[0].whr[0].im - Xi2[0].whr[0].im );
	 psi[site_c][3][jord].whr[0].im -= .5*( Xi1[0].whr[0].re - Xi2[0].whr[0].re );
	 psi[site_c][3][jord].whr[1].re += .5*( Xi1[0].whr[1].im - Xi2[0].whr[1].im );
	 psi[site_c][3][jord].whr[1].im -= .5*( Xi1[0].whr[1].re - Xi2[0].whr[1].re );
	 psi[site_c][3][jord].whr[2].re += .5*( Xi1[0].whr[2].im - Xi2[0].whr[2].im );
	 psi[site_c][3][jord].whr[2].im -= .5*( Xi1[0].whr[2].re - Xi2[0].whr[2].re );
	 
	 // mu = 1
	 point_up = Umu.Z->get(curr, 1, 1);
	 point_dn = Umu.Z->get(curr,-1, 1);
	 
	 Xi1[0].whr[0].re = psi[point_dn][0][kord].whr[0].re - psi[point_dn][3][kord].whr[0].re;    Xi1[0].whr[0].im = psi[point_dn][0][kord].whr[0].im - psi[point_dn][3][kord].whr[0].im;
	 Xi1[0].whr[1].re = psi[point_dn][0][kord].whr[1].re - psi[point_dn][3][kord].whr[1].re;    Xi1[0].whr[1].im = psi[point_dn][0][kord].whr[1].im - psi[point_dn][3][kord].whr[1].im;
	 Xi1[0].whr[2].re = psi[point_dn][0][kord].whr[2].re - psi[point_dn][3][kord].whr[2].re;    Xi1[0].whr[2].im = psi[point_dn][0][kord].whr[2].im - psi[point_dn][3][kord].whr[2].im;
	 Xi1[1].whr[0].re = psi[point_dn][1][kord].whr[0].re + psi[point_dn][2][kord].whr[0].re;    Xi1[1].whr[0].im = psi[point_dn][1][kord].whr[0].im + psi[point_dn][2][kord].whr[0].im;
	 Xi1[1].whr[1].re = psi[point_dn][1][kord].whr[1].re + psi[point_dn][2][kord].whr[1].re;    Xi1[1].whr[1].im = psi[point_dn][1][kord].whr[1].im + psi[point_dn][2][kord].whr[1].im;
	 Xi1[1].whr[2].re = psi[point_dn][1][kord].whr[2].re + psi[point_dn][2][kord].whr[2].re;    Xi1[1].whr[2].im = psi[point_dn][1][kord].whr[2].im + psi[point_dn][2][kord].whr[2].im;

	 Xi1[0] = dag(Umu.W[point_dn][1][jord-kord-1])*Xi1[0];
	 Xi1[1] = dag(Umu.W[point_dn][1][jord-kord-1])*Xi1[1];

	 Xi2[0].whr[0].re = psi[point_up][0][kord].whr[0].re + psi[point_up][3][kord].whr[0].re;    Xi2[0].whr[0].im = psi[point_up][0][kord].whr[0].im + psi[point_up][3][kord].whr[0].im;
	 Xi2[0].whr[1].re = psi[point_up][0][kord].whr[1].re + psi[point_up][3][kord].whr[1].re;    Xi2[0].whr[1].im = psi[point_up][0][kord].whr[1].im + psi[point_up][3][kord].whr[1].im;
	 Xi2[0].whr[2].re = psi[point_up][0][kord].whr[2].re + psi[point_up][3][kord].whr[2].re;    Xi2[0].whr[2].im = psi[point_up][0][kord].whr[2].im + psi[point_up][3][kord].whr[2].im;
	 Xi2[1].whr[0].re = psi[point_up][1][kord].whr[0].re - psi[point_up][2][kord].whr[0].re;    Xi2[1].whr[0].im = psi[point_up][1][kord].whr[0].im - psi[point_up][2][kord].whr[0].im;
	 Xi2[1].whr[1].re = psi[point_up][1][kord].whr[1].re - psi[point_up][2][kord].whr[1].re;    Xi2[1].whr[1].im = psi[point_up][1][kord].whr[1].im - psi[point_up][2][kord].whr[1].im;
	 Xi2[1].whr[2].re = psi[point_up][1][kord].whr[2].re - psi[point_up][2][kord].whr[2].re;    Xi2[1].whr[2].im = psi[point_up][1][kord].whr[2].im - psi[point_up][2][kord].whr[2].im;

	 Xi2[0] = Umu.W[curr][1][jord-kord-1]*Xi2[0];
	 Xi2[1] = Umu.W[curr][1][jord-kord-1]*Xi2[1];

	 // reconstruct
	 psi[site_c][0][jord] -= .5*(Xi1[0] + Xi2[0]);
	 psi[site_c][1][jord] -= .5*(Xi1[1] + Xi2[1]);

	 psi[site_c][2][jord] -= .5*( Xi1[1] - Xi2[1] );
	 psi[site_c][3][jord] += .5*( Xi1[0] - Xi2[0] );

	 
	 // mu = 2
	 point_up = Umu.Z->get(curr, 1, 2);
	 point_dn = Umu.Z->get(curr,-1, 2);

	 Xi1[0].whr[0].re = psi[point_dn][0][kord].whr[0].re + psi[point_dn][2][kord].whr[0].im;
	 Xi1[0].whr[0].im = psi[point_dn][0][kord].whr[0].im - psi[point_dn][2][kord].whr[0].re;
	 
	 Xi1[0].whr[1].re = psi[point_dn][0][kord].whr[1].re + psi[point_dn][2][kord].whr[1].im;
	 Xi1[0].whr[1].im = psi[point_dn][0][kord].whr[1].im - psi[point_dn][2][kord].whr[1].re;
	 
	 Xi1[0].whr[2].re = psi[point_dn][0][kord].whr[2].re + psi[point_dn][2][kord].whr[2].im;
	 Xi1[0].whr[2].im = psi[point_dn][0][kord].whr[2].im - psi[point_dn][2][kord].whr[2].re;
	 

	 Xi1[1].whr[0].re = psi[point_dn][1][kord].whr[0].re - psi[point_dn][3][kord].whr[0].im;
	 Xi1[1].whr[0].im = psi[point_dn][1][kord].whr[0].im + psi[point_dn][3][kord].whr[0].re;

	 Xi1[1].whr[1].re = psi[point_dn][1][kord].whr[1].re - psi[point_dn][3][kord].whr[1].im;
	 Xi1[1].whr[1].im = psi[point_dn][1][kord].whr[1].im + psi[point_dn][3][kord].whr[1].re;

	 Xi1[1].whr[2].re = psi[point_dn][1][kord].whr[2].re - psi[point_dn][3][kord].whr[2].im;
	 Xi1[1].whr[2].im = psi[point_dn][1][kord].whr[2].im + psi[point_dn][3][kord].whr[2].re;

	 Xi1[0] = dag(Umu.W[point_dn][2][jord-kord-1])*Xi1[0];
	 Xi1[1] = dag(Umu.W[point_dn][2][jord-kord-1])*Xi1[1];

	 Xi2[0].whr[0].re = psi[point_up][0][kord].whr[0].re - psi[point_up][2][kord].whr[0].im;
	 Xi2[0].whr[0].im = psi[point_up][0][kord].whr[0].im + psi[point_up][2][kord].whr[0].re;
	 
	 Xi2[0].whr[1].re = psi[point_up][0][kord].whr[1].re - psi[point_up][2][kord].whr[1].im;
	 Xi2[0].whr[1].im = psi[point_up][0][kord].whr[1].im + psi[point_up][2][kord].whr[1].re;
	 
	 Xi2[0].whr[2].re = psi[point_up][0][kord].whr[2].re - psi[point_up][2][kord].whr[2].im;
	 Xi2[0].whr[2].im = psi[point_up][0][kord].whr[2].im + psi[point_up][2][kord].whr[2].re;
	 
	 Xi2[1].whr[0].re = psi[point_up][1][kord].whr[0].re + psi[point_up][3][kord].whr[0].im;
	 Xi2[1].whr[0].im = psi[point_up][1][kord].whr[0].im - psi[point_up][3][kord].whr[0].re;
	 
	 Xi2[1].whr[1].re = psi[point_up][1][kord].whr[1].re + psi[point_up][3][kord].whr[1].im;
	 Xi2[1].whr[1].im = psi[point_up][1][kord].whr[1].im - psi[point_up][3][kord].whr[1].re;
	 
	 Xi2[1].whr[2].re = psi[point_up][1][kord].whr[2].re + psi[point_up][3][kord].whr[2].im;
	 Xi2[1].whr[2].im = psi[point_up][1][kord].whr[2].im - psi[point_up][3][kord].whr[2].re;

	 Xi2[0] = Umu.W[curr][2][jord-kord-1]*Xi2[0];
	 Xi2[1] = Umu.W[curr][2][jord-kord-1]*Xi2[1];

	 // reconstruct
	 psi[site_c][0][jord] -= .5*( Xi1[0] + Xi2[0] );
	 psi[site_c][1][jord] -= .5*( Xi1[1] + Xi2[1] );

	 psi[site_c][2][jord].whr[0].re += .5*( Xi1[0].whr[0].im - Xi2[0].whr[0].im );
	 psi[site_c][2][jord].whr[0].im -= .5*( Xi1[0].whr[0].re - Xi2[0].whr[0].re );
	 psi[site_c][2][jord].whr[1].re += .5*( Xi1[0].whr[1].im - Xi2[0].whr[1].im );
	 psi[site_c][2][jord].whr[1].im -= .5*( Xi1[0].whr[1].re - Xi2[0].whr[1].re );
	 psi[site_c][2][jord].whr[2].re += .5*( Xi1[0].whr[2].im - Xi2[0].whr[2].im );
	 psi[site_c][2][jord].whr[2].im -= .5*( Xi1[0].whr[2].re - Xi2[0].whr[2].re );

	 psi[site_c][3][jord].whr[0].re -= .5*( Xi1[1].whr[0].im - Xi2[1].whr[0].im );
	 psi[site_c][3][jord].whr[0].im += .5*( Xi1[1].whr[0].re - Xi2[1].whr[0].re );
	 psi[site_c][3][jord].whr[1].re -= .5*( Xi1[1].whr[1].im - Xi2[1].whr[1].im );
	 psi[site_c][3][jord].whr[1].im += .5*( Xi1[1].whr[1].re - Xi2[1].whr[1].re );
	 psi[site_c][3][jord].whr[2].re -= .5*( Xi1[1].whr[2].im - Xi2[1].whr[2].im );
	 psi[site_c][3][jord].whr[2].im += .5*( Xi1[1].whr[2].re - Xi2[1].whr[2].re );

	 // mu = 3
	 point_up = Umu.Z->get(curr, 1, 3);
	 point_dn = Umu.Z->get(curr,-1, 3);

	 Xi1[0].whr[0].re = psi[point_dn][0][kord].whr[0].re + psi[point_dn][2][kord].whr[0].re;    Xi1[0].whr[0].im = psi[point_dn][0][kord].whr[0].im + psi[point_dn][2][kord].whr[0].im;
	 Xi1[0].whr[1].re = psi[point_dn][0][kord].whr[1].re + psi[point_dn][2][kord].whr[1].re;    Xi1[0].whr[1].im = psi[point_dn][0][kord].whr[1].im + psi[point_dn][2][kord].whr[1].im;
	 Xi1[0].whr[2].re = psi[point_dn][0][kord].whr[2].re + psi[point_dn][2][kord].whr[2].re;    Xi1[0].whr[2].im = psi[point_dn][0][kord].whr[2].im + psi[point_dn][2][kord].whr[2].im;
	 Xi1[1].whr[0].re = psi[point_dn][1][kord].whr[0].re + psi[point_dn][3][kord].whr[0].re;    Xi1[1].whr[0].im = psi[point_dn][1][kord].whr[0].im + psi[point_dn][3][kord].whr[0].im;
	 Xi1[1].whr[1].re = psi[point_dn][1][kord].whr[1].re + psi[point_dn][3][kord].whr[1].re;    Xi1[1].whr[1].im = psi[point_dn][1][kord].whr[1].im + psi[point_dn][3][kord].whr[1].im;
	 Xi1[1].whr[2].re = psi[point_dn][1][kord].whr[2].re + psi[point_dn][3][kord].whr[2].re;    Xi1[1].whr[2].im = psi[point_dn][1][kord].whr[2].im + psi[point_dn][3][kord].whr[2].im;

	 Xi1[0] = dag(Umu.W[point_dn][3][jord-kord-1])*Xi1[0];
	 Xi1[1] = dag(Umu.W[point_dn][3][jord-kord-1])*Xi1[1];

	 Xi2[0].whr[0].re = psi[point_up][0][kord].whr[0].re - psi[point_up][2][kord].whr[0].re;    Xi2[0].whr[0].im = psi[point_up][0][kord].whr[0].im - psi[point_up][2][kord].whr[0].im;
	 Xi2[0].whr[1].re = psi[point_up][0][kord].whr[1].re - psi[point_up][2][kord].whr[1].re;    Xi2[0].whr[1].im = psi[point_up][0][kord].whr[1].im - psi[point_up][2][kord].whr[1].im;
	 Xi2[0].whr[2].re = psi[point_up][0][kord].whr[2].re - psi[point_up][2][kord].whr[2].re;    Xi2[0].whr[2].im = psi[point_up][0][kord].whr[2].im - psi[point_up][2][kord].whr[2].im;
	 Xi2[1].whr[0].re = psi[point_up][1][kord].whr[0].re - psi[point_up][3][kord].whr[0].re;    Xi2[1].whr[0].im = psi[point_up][1][kord].whr[0].im - psi[point_up][3][kord].whr[0].im;
	 Xi2[1].whr[1].re = psi[point_up][1][kord].whr[1].re - psi[point_up][3][kord].whr[1].re;    Xi2[1].whr[1].im = psi[point_up][1][kord].whr[1].im - psi[point_up][3][kord].whr[1].im;
	 Xi2[1].whr[2].re = psi[point_up][1][kord].whr[2].re - psi[point_up][3][kord].whr[2].re;    Xi2[1].whr[2].im = psi[point_up][1][kord].whr[2].im - psi[point_up][3][kord].whr[2].im;
	 
	 Xi2[0] = Umu.W[curr][3][jord-kord-1]*Xi2[0];
	 Xi2[1] = Umu.W[curr][3][jord-kord-1]*Xi2[1];
	 
	 // reconstruct
	 psi[site_c][0][jord] -= .5*(Xi1[0] + Xi2[0]);
	 psi[site_c][1][jord] -= .5*(Xi1[1] + Xi2[1]);
	 
	 psi[site_c][2][jord] -= .5*( Xi1[0] - Xi2[0] );
	 psi[site_c][3][jord] -= .5*( Xi1[1] - Xi2[1] );

#endif	  
	  
       } //kord
	
     } // curr
      
      
     // metto su scfld per applicare M0^-1, inverto e rimetto in psi
     for(int i = 0; i < Z->Size; i++){
       for(int mu = 0; mu < dim; mu++){
	 scfld->psi[i][mu] = psi[i][mu][jord];
       }
     }
      
     fft(0);
     M0inv();
     fft(1);
      
     for(int i = 0; i < Z->Size; i++){
       for(int mu = 0; mu < dim; mu++){
	 psi[i][mu][jord] = -(scfld->psi[i][mu]);
       }
     }
      
      
   }// jord

 }

#else

  // OpenMP Parallel

  void fillPT(int ptmax, ptGluon_fld &Umu){

    fft(0);
    M0inv();
    fft(1);

    int site_c, point_up, point_dn;
    SpinColor Xi1, Xi2;
    int curr,ttid;
    int chunk = (Z->Size)/NTHR;

    // per copiare non mi interessa l'ordine in cui accedo ai siti
#pragma omp parallel private(ttid) num_threads(NTHR)
    {
      ttid = omp_get_thread_num();
      for(int i = ttid*chunk; i < (ttid+1)*chunk; i++){
	for(int mu = 0; mu < dim; mu++){
	  psi[i][mu][0] = scfld->psi[i][mu];
	}
      } // siti
#pragma omp barrier
    } //parallel

    
    for(int jord = 1; jord <= ptmax; jord++){
      
#pragma omp parallel private(curr,site_c,ttid,point_up,point_dn,Xi1,Xi2) num_threads(NTHR)
      {
	ttid = omp_get_thread_num();
	
	for(int y0 = Z->Sz[0]*ttid/NTHR; y0 < Z->Sz[0]*(ttid+1)/NTHR; y0++){
	  for(int y1 = 0; y1 < Z->Sz[1]; y1++){
	    for(int y2 = 0; y2 < Z->Sz[2]; y2++){	  	      
	      for(int y3 = 0; y3 < Z->Sz[3]; y3++){

		curr = y3 + Z->Sz[3]*(y2 + Z->Sz[2]*(y1 + Z->Sz[1]*y0) );
		site_c = Z->get(curr);
		
		// moltiplicazione perturbativa per M
		for( int kord = 0; kord < jord; kord++) {
		  // massa critica
		  for (int mu = 0; mu < dim; mu++ ) {		  
		    psi[site_c][mu][jord] += mcpt[jord-kord]*psi[site_c][mu][kord];
		  }

#ifndef _USE_HALFSPINOR_
		  
		  for( int mu = 0; mu < dim; mu++){
		    point_up = Umu.Z->get(curr, 1, mu);
		    point_dn = Umu.Z->get(curr,-1, mu);
		    // 1 + \gamma_\mu

		    psi[point_dn].uno_p_gmu(Xi1,mu,kord);
		    psi[point_up].uno_m_gmu(Xi2,mu,kord);
		    
		    for( int nu = 0; nu < dim; nu++ ){
		      psi[site_c][nu][jord] -= .5*(dag(Umu.W[point_dn][mu][jord-kord-1])*Xi1[nu]);
		      psi[site_c][nu][jord] -= .5*( Umu.W[site_c][mu][jord-kord-1] )*Xi2[nu];
		    } //nu
		    
		  } // mu

#else	    
	 // -- Half Spinor --
	 // see arXiv:0905.3331v1 [hep-lat]

	 // mu = 0
	 point_up = Umu.Z->get(curr, 1, 0);
	 point_dn = Umu.Z->get(curr,-1, 0);
	 
	 Xi1[0].whr[0].re = psi[point_dn][0][kord].whr[0].re + psi[point_dn][3][kord].whr[0].im;    Xi1[0].whr[0].im = psi[point_dn][0][kord].whr[0].im - psi[point_dn][3][kord].whr[0].re;
	 Xi1[0].whr[1].re = psi[point_dn][0][kord].whr[1].re + psi[point_dn][3][kord].whr[1].im;    Xi1[0].whr[1].im = psi[point_dn][0][kord].whr[1].im - psi[point_dn][3][kord].whr[1].re;
	 Xi1[0].whr[2].re = psi[point_dn][0][kord].whr[2].re + psi[point_dn][3][kord].whr[2].im;    Xi1[0].whr[2].im = psi[point_dn][0][kord].whr[2].im - psi[point_dn][3][kord].whr[2].re;
	 Xi1[1].whr[0].re = psi[point_dn][1][kord].whr[0].re + psi[point_dn][2][kord].whr[0].im;    Xi1[1].whr[0].im = psi[point_dn][1][kord].whr[0].im - psi[point_dn][2][kord].whr[0].re;
	 Xi1[1].whr[1].re = psi[point_dn][1][kord].whr[1].re + psi[point_dn][2][kord].whr[1].im;    Xi1[1].whr[1].im = psi[point_dn][1][kord].whr[1].im - psi[point_dn][2][kord].whr[1].re;
	 Xi1[1].whr[2].re = psi[point_dn][1][kord].whr[2].re + psi[point_dn][2][kord].whr[2].im;    Xi1[1].whr[2].im = psi[point_dn][1][kord].whr[2].im - psi[point_dn][2][kord].whr[2].re;

	 Xi1[0] = dag(Umu.W[point_dn][0][jord-kord-1])*Xi1[0];
	 Xi1[1] = dag(Umu.W[point_dn][0][jord-kord-1])*Xi1[1];

	 Xi2[0].whr[0].re = psi[point_up][0][kord].whr[0].re - psi[point_up][3][kord].whr[0].im;      Xi2[0].whr[0].im = psi[point_up][0][kord].whr[0].im + psi[point_up][3][kord].whr[0].re;
	 Xi2[0].whr[1].re = psi[point_up][0][kord].whr[1].re - psi[point_up][3][kord].whr[1].im;      Xi2[0].whr[1].im = psi[point_up][0][kord].whr[1].im + psi[point_up][3][kord].whr[1].re;
	 Xi2[0].whr[2].re = psi[point_up][0][kord].whr[2].re - psi[point_up][3][kord].whr[2].im;      Xi2[0].whr[2].im = psi[point_up][0][kord].whr[2].im + psi[point_up][3][kord].whr[2].re;
	 Xi2[1].whr[0].re = psi[point_up][1][kord].whr[0].re - psi[point_up][2][kord].whr[0].im;      Xi2[1].whr[0].im = psi[point_up][1][kord].whr[0].im + psi[point_up][2][kord].whr[0].re;
	 Xi2[1].whr[1].re = psi[point_up][1][kord].whr[1].re - psi[point_up][2][kord].whr[1].im;      Xi2[1].whr[1].im = psi[point_up][1][kord].whr[1].im + psi[point_up][2][kord].whr[1].re;
	 Xi2[1].whr[2].re = psi[point_up][1][kord].whr[2].re - psi[point_up][2][kord].whr[2].im;      Xi2[1].whr[2].im = psi[point_up][1][kord].whr[2].im + psi[point_up][2][kord].whr[2].re;

	 Xi2[0] = Umu.W[curr][0][jord-kord-1]*Xi2[0];
	 Xi2[1] = Umu.W[curr][0][jord-kord-1]*Xi2[1];
	 
	 // reconstruct
	 psi[site_c][0][jord] -= .5*(Xi1[0] + Xi2[0]);
	 psi[site_c][1][jord] -= .5*(Xi1[1] + Xi2[1]);

	 psi[site_c][2][jord].whr[0].re += .5*( Xi1[1].whr[0].im - Xi2[1].whr[0].im );
	 psi[site_c][2][jord].whr[0].im -= .5*( Xi1[1].whr[0].re - Xi2[1].whr[0].re );
	 psi[site_c][2][jord].whr[1].re += .5*( Xi1[1].whr[1].im - Xi2[1].whr[1].im );
	 psi[site_c][2][jord].whr[1].im -= .5*( Xi1[1].whr[1].re - Xi2[1].whr[1].re );
	 psi[site_c][2][jord].whr[2].re += .5*( Xi1[1].whr[2].im - Xi2[1].whr[2].im );
	 psi[site_c][2][jord].whr[2].im -= .5*( Xi1[1].whr[2].re - Xi2[1].whr[2].re );

	 psi[site_c][3][jord].whr[0].re += .5*( Xi1[0].whr[0].im - Xi2[0].whr[0].im );
	 psi[site_c][3][jord].whr[0].im -= .5*( Xi1[0].whr[0].re - Xi2[0].whr[0].re );
	 psi[site_c][3][jord].whr[1].re += .5*( Xi1[0].whr[1].im - Xi2[0].whr[1].im );
	 psi[site_c][3][jord].whr[1].im -= .5*( Xi1[0].whr[1].re - Xi2[0].whr[1].re );
	 psi[site_c][3][jord].whr[2].re += .5*( Xi1[0].whr[2].im - Xi2[0].whr[2].im );
	 psi[site_c][3][jord].whr[2].im -= .5*( Xi1[0].whr[2].re - Xi2[0].whr[2].re );
	 
	 // mu = 1
	 point_up = Umu.Z->get(curr, 1, 1);
	 point_dn = Umu.Z->get(curr,-1, 1);
	 
	 Xi1[0].whr[0].re = psi[point_dn][0][kord].whr[0].re - psi[point_dn][3][kord].whr[0].re;    Xi1[0].whr[0].im = psi[point_dn][0][kord].whr[0].im - psi[point_dn][3][kord].whr[0].im;
	 Xi1[0].whr[1].re = psi[point_dn][0][kord].whr[1].re - psi[point_dn][3][kord].whr[1].re;    Xi1[0].whr[1].im = psi[point_dn][0][kord].whr[1].im - psi[point_dn][3][kord].whr[1].im;
	 Xi1[0].whr[2].re = psi[point_dn][0][kord].whr[2].re - psi[point_dn][3][kord].whr[2].re;    Xi1[0].whr[2].im = psi[point_dn][0][kord].whr[2].im - psi[point_dn][3][kord].whr[2].im;
	 Xi1[1].whr[0].re = psi[point_dn][1][kord].whr[0].re + psi[point_dn][2][kord].whr[0].re;    Xi1[1].whr[0].im = psi[point_dn][1][kord].whr[0].im + psi[point_dn][2][kord].whr[0].im;
	 Xi1[1].whr[1].re = psi[point_dn][1][kord].whr[1].re + psi[point_dn][2][kord].whr[1].re;    Xi1[1].whr[1].im = psi[point_dn][1][kord].whr[1].im + psi[point_dn][2][kord].whr[1].im;
	 Xi1[1].whr[2].re = psi[point_dn][1][kord].whr[2].re + psi[point_dn][2][kord].whr[2].re;    Xi1[1].whr[2].im = psi[point_dn][1][kord].whr[2].im + psi[point_dn][2][kord].whr[2].im;

	 Xi1[0] = dag(Umu.W[point_dn][1][jord-kord-1])*Xi1[0];
	 Xi1[1] = dag(Umu.W[point_dn][1][jord-kord-1])*Xi1[1];

	 Xi2[0].whr[0].re = psi[point_up][0][kord].whr[0].re + psi[point_up][3][kord].whr[0].re;    Xi2[0].whr[0].im = psi[point_up][0][kord].whr[0].im + psi[point_up][3][kord].whr[0].im;
	 Xi2[0].whr[1].re = psi[point_up][0][kord].whr[1].re + psi[point_up][3][kord].whr[1].re;    Xi2[0].whr[1].im = psi[point_up][0][kord].whr[1].im + psi[point_up][3][kord].whr[1].im;
	 Xi2[0].whr[2].re = psi[point_up][0][kord].whr[2].re + psi[point_up][3][kord].whr[2].re;    Xi2[0].whr[2].im = psi[point_up][0][kord].whr[2].im + psi[point_up][3][kord].whr[2].im;
	 Xi2[1].whr[0].re = psi[point_up][1][kord].whr[0].re - psi[point_up][2][kord].whr[0].re;    Xi2[1].whr[0].im = psi[point_up][1][kord].whr[0].im - psi[point_up][2][kord].whr[0].im;
	 Xi2[1].whr[1].re = psi[point_up][1][kord].whr[1].re - psi[point_up][2][kord].whr[1].re;    Xi2[1].whr[1].im = psi[point_up][1][kord].whr[1].im - psi[point_up][2][kord].whr[1].im;
	 Xi2[1].whr[2].re = psi[point_up][1][kord].whr[2].re - psi[point_up][2][kord].whr[2].re;    Xi2[1].whr[2].im = psi[point_up][1][kord].whr[2].im - psi[point_up][2][kord].whr[2].im;

	 Xi2[0] = Umu.W[curr][1][jord-kord-1]*Xi2[0];
	 Xi2[1] = Umu.W[curr][1][jord-kord-1]*Xi2[1];

	 // reconstruct
	 psi[site_c][0][jord] -= .5*(Xi1[0] + Xi2[0]);
	 psi[site_c][1][jord] -= .5*(Xi1[1] + Xi2[1]);

	 psi[site_c][2][jord] -= .5*( Xi1[1] - Xi2[1] );
	 psi[site_c][3][jord] += .5*( Xi1[0] - Xi2[0] );

	 
	 // mu = 2
	 point_up = Umu.Z->get(curr, 1, 2);
	 point_dn = Umu.Z->get(curr,-1, 2);

	 Xi1[0].whr[0].re = psi[point_dn][0][kord].whr[0].re + psi[point_dn][2][kord].whr[0].im;
	 Xi1[0].whr[0].im = psi[point_dn][0][kord].whr[0].im - psi[point_dn][2][kord].whr[0].re;
	 
	 Xi1[0].whr[1].re = psi[point_dn][0][kord].whr[1].re + psi[point_dn][2][kord].whr[1].im;
	 Xi1[0].whr[1].im = psi[point_dn][0][kord].whr[1].im - psi[point_dn][2][kord].whr[1].re;
	 
	 Xi1[0].whr[2].re = psi[point_dn][0][kord].whr[2].re + psi[point_dn][2][kord].whr[2].im;
	 Xi1[0].whr[2].im = psi[point_dn][0][kord].whr[2].im - psi[point_dn][2][kord].whr[2].re;
	 

	 Xi1[1].whr[0].re = psi[point_dn][1][kord].whr[0].re - psi[point_dn][3][kord].whr[0].im;
	 Xi1[1].whr[0].im = psi[point_dn][1][kord].whr[0].im + psi[point_dn][3][kord].whr[0].re;

	 Xi1[1].whr[1].re = psi[point_dn][1][kord].whr[1].re - psi[point_dn][3][kord].whr[1].im;
	 Xi1[1].whr[1].im = psi[point_dn][1][kord].whr[1].im + psi[point_dn][3][kord].whr[1].re;

	 Xi1[1].whr[2].re = psi[point_dn][1][kord].whr[2].re - psi[point_dn][3][kord].whr[2].im;
	 Xi1[1].whr[2].im = psi[point_dn][1][kord].whr[2].im + psi[point_dn][3][kord].whr[2].re;

	 Xi1[0] = dag(Umu.W[point_dn][2][jord-kord-1])*Xi1[0];
	 Xi1[1] = dag(Umu.W[point_dn][2][jord-kord-1])*Xi1[1];

	 Xi2[0].whr[0].re = psi[point_up][0][kord].whr[0].re - psi[point_up][2][kord].whr[0].im;
	 Xi2[0].whr[0].im = psi[point_up][0][kord].whr[0].im + psi[point_up][2][kord].whr[0].re;
	 
	 Xi2[0].whr[1].re = psi[point_up][0][kord].whr[1].re - psi[point_up][2][kord].whr[1].im;
	 Xi2[0].whr[1].im = psi[point_up][0][kord].whr[1].im + psi[point_up][2][kord].whr[1].re;
	 
	 Xi2[0].whr[2].re = psi[point_up][0][kord].whr[2].re - psi[point_up][2][kord].whr[2].im;
	 Xi2[0].whr[2].im = psi[point_up][0][kord].whr[2].im + psi[point_up][2][kord].whr[2].re;
	 
	 Xi2[1].whr[0].re = psi[point_up][1][kord].whr[0].re + psi[point_up][3][kord].whr[0].im;
	 Xi2[1].whr[0].im = psi[point_up][1][kord].whr[0].im - psi[point_up][3][kord].whr[0].re;
	 
	 Xi2[1].whr[1].re = psi[point_up][1][kord].whr[1].re + psi[point_up][3][kord].whr[1].im;
	 Xi2[1].whr[1].im = psi[point_up][1][kord].whr[1].im - psi[point_up][3][kord].whr[1].re;
	 
	 Xi2[1].whr[2].re = psi[point_up][1][kord].whr[2].re + psi[point_up][3][kord].whr[2].im;
	 Xi2[1].whr[2].im = psi[point_up][1][kord].whr[2].im - psi[point_up][3][kord].whr[2].re;

	 Xi2[0] = Umu.W[curr][2][jord-kord-1]*Xi2[0];
	 Xi2[1] = Umu.W[curr][2][jord-kord-1]*Xi2[1];

	 // reconstruct
	 psi[site_c][0][jord] -= .5*( Xi1[0] + Xi2[0] );
	 psi[site_c][1][jord] -= .5*( Xi1[1] + Xi2[1] );

	 psi[site_c][2][jord].whr[0].re += .5*( Xi1[0].whr[0].im - Xi2[0].whr[0].im );
	 psi[site_c][2][jord].whr[0].im -= .5*( Xi1[0].whr[0].re - Xi2[0].whr[0].re );
	 psi[site_c][2][jord].whr[1].re += .5*( Xi1[0].whr[1].im - Xi2[0].whr[1].im );
	 psi[site_c][2][jord].whr[1].im -= .5*( Xi1[0].whr[1].re - Xi2[0].whr[1].re );
	 psi[site_c][2][jord].whr[2].re += .5*( Xi1[0].whr[2].im - Xi2[0].whr[2].im );
	 psi[site_c][2][jord].whr[2].im -= .5*( Xi1[0].whr[2].re - Xi2[0].whr[2].re );

	 psi[site_c][3][jord].whr[0].re -= .5*( Xi1[1].whr[0].im - Xi2[1].whr[0].im );
	 psi[site_c][3][jord].whr[0].im += .5*( Xi1[1].whr[0].re - Xi2[1].whr[0].re );
	 psi[site_c][3][jord].whr[1].re -= .5*( Xi1[1].whr[1].im - Xi2[1].whr[1].im );
	 psi[site_c][3][jord].whr[1].im += .5*( Xi1[1].whr[1].re - Xi2[1].whr[1].re );
	 psi[site_c][3][jord].whr[2].re -= .5*( Xi1[1].whr[2].im - Xi2[1].whr[2].im );
	 psi[site_c][3][jord].whr[2].im += .5*( Xi1[1].whr[2].re - Xi2[1].whr[2].re );

	 // mu = 3
	 point_up = Umu.Z->get(curr, 1, 3);
	 point_dn = Umu.Z->get(curr,-1, 3);

	 Xi1[0].whr[0].re = psi[point_dn][0][kord].whr[0].re + psi[point_dn][2][kord].whr[0].re;    Xi1[0].whr[0].im = psi[point_dn][0][kord].whr[0].im + psi[point_dn][2][kord].whr[0].im;
	 Xi1[0].whr[1].re = psi[point_dn][0][kord].whr[1].re + psi[point_dn][2][kord].whr[1].re;    Xi1[0].whr[1].im = psi[point_dn][0][kord].whr[1].im + psi[point_dn][2][kord].whr[1].im;
	 Xi1[0].whr[2].re = psi[point_dn][0][kord].whr[2].re + psi[point_dn][2][kord].whr[2].re;    Xi1[0].whr[2].im = psi[point_dn][0][kord].whr[2].im + psi[point_dn][2][kord].whr[2].im;
	 Xi1[1].whr[0].re = psi[point_dn][1][kord].whr[0].re + psi[point_dn][3][kord].whr[0].re;    Xi1[1].whr[0].im = psi[point_dn][1][kord].whr[0].im + psi[point_dn][3][kord].whr[0].im;
	 Xi1[1].whr[1].re = psi[point_dn][1][kord].whr[1].re + psi[point_dn][3][kord].whr[1].re;    Xi1[1].whr[1].im = psi[point_dn][1][kord].whr[1].im + psi[point_dn][3][kord].whr[1].im;
	 Xi1[1].whr[2].re = psi[point_dn][1][kord].whr[2].re + psi[point_dn][3][kord].whr[2].re;    Xi1[1].whr[2].im = psi[point_dn][1][kord].whr[2].im + psi[point_dn][3][kord].whr[2].im;

	 Xi1[0] = dag(Umu.W[point_dn][3][jord-kord-1])*Xi1[0];
	 Xi1[1] = dag(Umu.W[point_dn][3][jord-kord-1])*Xi1[1];

	 Xi2[0].whr[0].re = psi[point_up][0][kord].whr[0].re - psi[point_up][2][kord].whr[0].re;    Xi2[0].whr[0].im = psi[point_up][0][kord].whr[0].im - psi[point_up][2][kord].whr[0].im;
	 Xi2[0].whr[1].re = psi[point_up][0][kord].whr[1].re - psi[point_up][2][kord].whr[1].re;    Xi2[0].whr[1].im = psi[point_up][0][kord].whr[1].im - psi[point_up][2][kord].whr[1].im;
	 Xi2[0].whr[2].re = psi[point_up][0][kord].whr[2].re - psi[point_up][2][kord].whr[2].re;    Xi2[0].whr[2].im = psi[point_up][0][kord].whr[2].im - psi[point_up][2][kord].whr[2].im;
	 Xi2[1].whr[0].re = psi[point_up][1][kord].whr[0].re - psi[point_up][3][kord].whr[0].re;    Xi2[1].whr[0].im = psi[point_up][1][kord].whr[0].im - psi[point_up][3][kord].whr[0].im;
	 Xi2[1].whr[1].re = psi[point_up][1][kord].whr[1].re - psi[point_up][3][kord].whr[1].re;    Xi2[1].whr[1].im = psi[point_up][1][kord].whr[1].im - psi[point_up][3][kord].whr[1].im;
	 Xi2[1].whr[2].re = psi[point_up][1][kord].whr[2].re - psi[point_up][3][kord].whr[2].re;    Xi2[1].whr[2].im = psi[point_up][1][kord].whr[2].im - psi[point_up][3][kord].whr[2].im;
	 
	 Xi2[0] = Umu.W[curr][3][jord-kord-1]*Xi2[0];
	 Xi2[1] = Umu.W[curr][3][jord-kord-1]*Xi2[1];
	 
	 // reconstruct
	 psi[site_c][0][jord] -= .5*(Xi1[0] + Xi2[0]);
	 psi[site_c][1][jord] -= .5*(Xi1[1] + Xi2[1]);
	 
	 psi[site_c][2][jord] -= .5*( Xi1[0] - Xi2[0] );
	 psi[site_c][3][jord] -= .5*( Xi1[1] - Xi2[1] );

#endif	  

		} //kord
	      } // siti
	    }
	  }
	}
#pragma omp barrier
      
      // metto su scfld per applicare M0^-1, inverto e rimetto in psi
      // per copiare non mi interessa l'ordine in cui accedo ai siti
	for(int i = ttid*chunk; i < (ttid+1)*chunk; i++){
	  for(int mu = 0; mu < dim; mu++){
	    scfld->psi[i][mu] = psi[i][mu][jord];
	  }
	} // siti
#pragma omp barrier
      } //parallel
 
      
      fft(0);
      M0inv();
      fft(1);
      
      // per copiare non mi interessa l'ordine in cui accedo ai siti
#pragma omp parallel private(ttid) num_threads(NTHR)
      {
	ttid = omp_get_thread_num();
	for(int i = ttid*chunk; i < (ttid+1)*chunk; i++){
	  for(int mu = 0; mu < dim; mu++){
	    psi[i][mu][jord] = -(scfld->psi[i][mu]);
	  }
	} // siti
#pragma omp barrier
      } //parallel
    }// jord

  }
#endif  


};


#endif
