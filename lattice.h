#include <stdlib.h>

#define dim 4
#define x0 1
#define x1 2
#define x2 3
#define x3 4

#define MAX(x, y) ((x) > (y) ? (x) : (y))
#define MAX4(x, y, z, t) (MAX(MAX((x), (y)), MAX((z), (t))))

#define GETINDEX(x, y, z, t) ((x)*(sizeYZT) + (y)*(sizeZT) + (z)*(sizeT) + (t))
#define MOD(i, m) (((m) + (i)) % (m))

#define PBC

const int lat_offset = dim+1;

class latt;

class SU3_fld;
class Gluon_fld;
class CVector_fld;
class SpinColor_fld;

class ptSU3_fld;
class ptGluon_fld;
class ptCVector_fld;
class ptSpinColor_fld;

// Important about the current structure of class lattice: the mail structure
// is the array of array LL. This is a vector of size Size, and eache element
// is an array of dimension 9. The order of the array is the "phisical" order
// of space-time. The 5th(#4) element of the array is the order in which the
// sites of the fields are stored in computer memory. The elements form 0 to 3
// are the locations of the previous next neighbours, while the elements from
// 5 to 8 are the successive. (See down)
//
//    #0    #1    #2    #3    #4    #5  #6   #7   #8
//  .     .                 .     .                   .
//  .......                 -------                   .
//  |     |                 | n-1 |                   |
//  |-----------------------|-----|--------------------
//  | -x0 | -x1 | -x2 | -x3 | n   | x0 | x1 | x2 | x3 |
//  |-----------------------|-----|--------------------
//  |                       | n+1 |                   |
//  .                       |-----|                   .
//  .                       .     .                   .
//


class latt{
  friend class SU3_fld;
  friend class Gluon_fld;
  friend int get(SU3_fld*, int *);  
  friend int get(Gluon_fld *, int *);

  friend class CVector_fld;
  friend class SpinColor_fld;
  friend int get(CVector_fld*, int *);  
  friend int get(SpinColor_fld *, int *);

  friend class ptSU3_fld;
  friend class ptGluon_fld;
  friend class ptCVector_fld;
  friend class ptSpinColor_fld;

  friend int get(ptSU3_fld*, int *);  
  friend int get(ptSU3_fld*, int, int, int);
  friend int get(ptSU3_fld*, int);

  friend int get(ptGluon_fld *, int *);
  friend ptSU3 staple(int n, int mu, int nu);
  friend int get(ptCVector_fld*, int *);  
  friend int get(ptCVector_fld*, int, int, int);
  friend int get(ptCVector_fld*, int);

  friend int get(ptSpinColor_fld *, int *);

public:
  int *Sz;

#if (dim == 4)

  int sizeYZT; 
  int sizeZT;
  int sizeT;
  int Size;   

  double **p;
  double **pbar;
  double **p2bar;
  double **p2hat;

#endif


  //int **L;  
  std::vector<std::vector<int> > L;
  latt(int *size) : L(size[0]*size[1]*size[2]*size[3], 
                      std::vector<int>(1 + 2*dim) ){
    int c = 0;

#if (dim == 4)
    Sz = new int[4];
    Sz[0] = size[0];
    Sz[1] = size[1];
    Sz[2] = size[2];
    Sz[3] = size[3];
    sizeT = Sz[3];
    sizeZT= Sz[2]*Sz[3];
    sizeYZT  = Sz[1]*Sz[2]*Sz[3];
    Size    = Sz[0]*Sz[1]*Sz[2]*Sz[3];
    
    //L = new int*[Size];
    //
    //for(int ii=0; ii<Size; ii++){
    //  L[ii] = new int[1+2*dim];      
    //}    

    for (int i = 0; i < Size; i++){
      // il sito nel reticolo deve corrispondere alla locazione fisica, mentre
      // il numero contenuto deve essere la locazione di memoria
      L[i][dim] = c++;
      //printf("%d  ->%d \n", i,L[i][4]);
    }  


    int i;
    for (int x = 0; x < Sz[0]; x++)
	for (int y = 0; y < Sz[1]; y++)
	    for (int z = 0; z < Sz[2]; z++)
		for (int t = 0; t < Sz[3]; t++) {
		    i = GETINDEX(x,y,z,t);

		    L[i][0] = GETINDEX(MOD(x-1,Sz[0]), y, z, t);
		    L[i][1] = GETINDEX(x, MOD(y-1,Sz[1]), z, t);
		    L[i][2] = GETINDEX(x, y, MOD(z-1,Sz[2]), t);
		    L[i][3] = GETINDEX(x, y, z, MOD(t-1,Sz[3]));

		    L[i][5] = GETINDEX((x+1)%Sz[0], y, z, t);
		    L[i][6] = GETINDEX(x, (y+1)%Sz[1], z, t);
		    L[i][7] = GETINDEX(x, y, (z+1)%Sz[2], t);
		    L[i][8] = GETINDEX(x, y, z, (t+1)%Sz[3]);
		    
		}
  }

  

  void p_init(){
    p     = new double* [4];
    pbar  = new double* [4];
    p2bar = new double* [4];
    p2hat = new double* [4];

    for (int i = 0; i < dim; i++) {
	p[i]     = new double [Sz[i]];
	pbar[i]  = new double [Sz[i]];
	p2bar[i] = new double [Sz[i]];
	p2hat[i] = new double [Sz[i]];
    }

#ifdef ABC
    for (int i = 0; i < Sz[0]; i++) {
      p[0][i] = 2 * M_PI * (i+.5) / Sz[0];
      pbar[0][i] = sin(p[0][i]);
      p2bar[0][i] = pbar[0][i]*pbar[0][i];
      p2hat[0][i] = sin(p[0][i] * 0.5); p2hat[0][i] *= 4*p2hat[0][i];
    }
#elif defined PBC
    for (int i = 0; i < Sz[0]; i++) {
      p[0][i] = 2 * M_PI * i / Sz[0];
      pbar[0][i] = sin(p[0][i]);
      p2bar[0][i] = pbar[0][i]*pbar[0][i];
      p2hat[0][i] = sin(p[0][i] * 0.5); p2hat[0][i] *= 4*p2hat[0][i];
    }
#endif
    for (int i = 0; i < Sz[1]; i++) {
	p[1][i] = 2 * M_PI * i / Sz[1];
	pbar[1][i] = sin(p[1][i]);
	p2bar[1][i] = pbar[1][i]*pbar[1][i];
	p2hat[1][i] = sin(p[1][i] * 0.5); p2hat[1][i] *= 4*p2hat[1][i];
    }
    for (int i = 0; i < Sz[2]; i++) {
	p[2][i] = 2 * M_PI * i / Sz[2];
	pbar[2][i] = sin(p[2][i]);
	p2bar[2][i] = pbar[2][i]*pbar[2][i];
	p2hat[2][i] = sin(p[2][i] * 0.5); p2hat[2][i] *= 4*p2hat[2][i];
    }
    for (int i = 0; i < Sz[3]; i++) {
	p[3][i] = 2 * M_PI * i / Sz[3];
	pbar[3][i] = sin(p[3][i]);
	p2bar[3][i] = pbar[3][i]*pbar[3][i];
	p2hat[3][i] = sin(p[3][i] * 0.5); p2hat[3][i] *= 4*p2hat[3][i];
    }

#endif
  }


  void p_pbc(){
    for (int i = 0; i < Sz[0]; i++) {
      p[0][i] = 2 * M_PI * i / Sz[0];
      pbar[0][i] = sin(p[0][i]);
      p2bar[0][i] = pbar[0][i]*pbar[0][i];
      p2hat[0][i] = sin(p[0][i] * 0.5); p2hat[0][i] *= 4*p2hat[0][i];
    }
  }

  void p_abc(){
    for (int i = 0; i < Sz[0]; i++) {
      p[0][i] = 2 * M_PI * (i+.5) / Sz[0];
      pbar[0][i] = sin(p[0][i]);
      p2bar[0][i] = pbar[0][i]*pbar[0][i];
      p2hat[0][i] = sin(p[0][i] * 0.5); p2hat[0][i] *= 4*p2hat[0][i];
    }
  }


  int get(int *n) {
      return (GETINDEX(n[0], n[1], n[2], n[3])); 
  }

  void get(int i, int *coord) {
      int modYZT, modZT;

      coord[0] = i / sizeYZT;
      coord[1] = (modYZT = i % sizeYZT) / sizeZT;
      coord[2] = (modZT = modYZT % sizeZT) / sizeT;
      coord[3] = modZT % sizeT;
  }
  
  int get(int n) {
    return L[n][dim];
  }
  
  int get(int n, int step, int dir) {
      if (step < 0)
	  for (int i = 0; i < (-1 * step); i++)
	      n = L[n][dir];
      else
	  for (int i = 0; i < step; i++)
	      n = L[n][lat_offset+dir];
      
      return L[n][dim];
  }



  int get(int n, int step1, int dir1, int step2, int dir2) {

    if(step1 < 0){
      for(int i = 0; i < abs(step1); i++){
	n = L[n][dir1];
      }
    }
    else{
      for(int i = 0; i < step1; i++){
	n = L[n][lat_offset+dir1];
      }
    }

    if(step2 < 0){
      for(int i = 0; i < abs(step2); i++){
	n = L[n][dir2];
      }
    }
    else{
      for(int i = 0; i < step2; i++){
	n = L[n][lat_offset+dir2];
      }
    }

    return L[n][dim];
  }

  int near(int n, int sign, int dir){
      if (sign < 0)
	  return L[n][dir]+1;
      else
	  return L[n][lat_offset+dir]+1;
  }

  int near2(int n, int sign, int dir){
      if ((dir == 0 && sign == +1 && (n / sizeYZT) == Sz[0]-1) || (dir == 0 && sign == -1 && (n / sizeYZT) == 0)) {
	  if (sign < 0)
	      return -(L[n][dir]+1);
	  else
	      return -(L[n][lat_offset+dir]+1);
      }
      else {
	  if (sign < 0)
	      return +(L[n][dir]+1);
	  else
	      return +(L[n][lat_offset+dir]+1);
      }
  }



};    // end class latt

