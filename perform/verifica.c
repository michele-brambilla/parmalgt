//---------------------------------------------|
//Include Section ---------------------- BEGIN |
//                                             |
//                                             |
#undef RANF
//  Standard libs                              |
#include <iostream>
#include<cmath>

//  algebra classes                            |
//#include "MyMath.h"
#include<QCDpt.h>
//                                             |
//                                             |
//Include Section ---------------------- END   |
//---------------------------------------------|



#include<ctime>
#include<sys/time.h>
#include<fstream>
#include<stdlib.h>


typedef struct{
  int i;
  double flops;
  double flops2;
  double tempo;
  double tempo2;
} t_data;

using namespace std;

int main(int argc, char *argv[]) {

  int START = 2;
  int STOP = 5;
  int RIPETIZIONI = 5;
  double oo;

  
  /* per clock() */
  t_data data[STOP -START +1];
  clock_t t1, t2;

  /* per gettimeoday */
  struct timeval tv1, tv2;
  //  unsigned   int o;

  ofstream fout;
  fout.open("./result/res.dat");
  cout.precision(5);

  /* Alloca lunghezza ripetizioni */
  for(int j = 0; j<=STOP-START; j++){
    data[j].i = (int)pow(10.,START+j);
    cout << data[j].i << endl;
  }


  SU3 A,B,C;
  /* Alloca matrici scalari */
  srand(time(NULL));
  for(int i = 0; i < 9; i++){
    A.whr[i].re = rand();
    A.whr[i].im = rand();
    B.whr[i].re = rand();
    B.whr[i].im = rand();
  }

  fout << "#iteraz\ttempo\tflops\n" << endl;
  fout << "#index = 0 - addizione matrici scalari #\n" << endl;

  /* Calcolo fps: addizioni matrici scalari */  
  oo = (2*(9*2));
  for(int j = 0; j <= STOP-START; j++){
    
    data[j].flops  = 0;
    data[j].flops2 = 0;
    data[j].tempo  = 0;
    data[j].tempo2 = 0;

    for(int r = 0; r < RIPETIZIONI; r++){
      gettimeofday(&tv1, NULL);
      t1 = clock();
      for(int i = 0; i < data[j].i ; i++){
	A += B;
	A -= B;
      }
      gettimeofday(&tv2, NULL);
      t2 = clock() - t1;

      printf("%f\n",A.whr[0].re);
      data[j].tempo   += t2/( double)CLOCKS_PER_SEC;
      data[j].tempo2  += ((tv2.tv_sec - tv1.tv_sec) + 
			  (tv2.tv_usec - tv1.tv_usec)/( double)1000000);
    } 
    
    data[j].flops   = oo*RIPETIZIONI*data[j].i / data[j].tempo ;
    data[j].flops2  = oo*RIPETIZIONI*data[j].i / data[j].tempo2;

    fout << START+j         << "\t" 
	 << scientific << data[j].tempo   << "\t" 
	 << scientific << data[j].tempo2  << "\t"
	 << scientific << data[j].flops*1e-9   << "\t" 
	 << scientific << data[j].flops2*1e-9  << endl;
  }
  

  /* Alloca matrici perturbative */
  ptSU3 AA,BB;
  for(int i = 0; i < PTORD; i++){
    AA.ptU[i] = A;
    BB.ptU[i] = B;
  }

  /* Calcolo fps: addizioni matrici perturbative */
  oo = (2*(2+PTORD*9*2));
  fout << "\n\n#index = 1 - addizione matrici perturbative #\n" << endl;
  for(int j = 0; j <= STOP-START; j++){

    data[j].flops  = 0;
    data[j].flops2 = 0;
    data[j].tempo  = 0;
    data[j].tempo2 = 0;

    for(int r = 0; r < RIPETIZIONI; r++){
      
      gettimeofday(&tv1, NULL);
      t1 = clock();
      for(int i = 0; i < data[j].i ; i++){
	AA += BB;
	AA -= BB;
      }
      gettimeofday(&tv2, NULL);
      t2 = clock() - t1;
      printf("%f\n",AA.ptU[2].whr[0].re);
      data[j].tempo   += t2/( double)CLOCKS_PER_SEC;
      data[j].tempo2  += ((tv2.tv_sec - tv1.tv_sec) +
			  (tv2.tv_usec - tv1.tv_usec)/( double)1000000);
    }
    data[j].flops   = oo*RIPETIZIONI*data[j].i / data[j].tempo;
    data[j].flops2  = oo*RIPETIZIONI*data[j].i / data[j].tempo2;

    fout << START+j         << "\t"
	 << scientific << data[j].tempo   << "\t"
	 << scientific << data[j].tempo2  << "\t"
	 << scientific << data[j].flops*1e-9   << "\t"
	 << scientific << data[j].flops2*1e-9  << endl;

  }


  /* Alloca matrici per moltiplicazioni */
  B.whr[0].re = 1;
  B.whr[1].re = 2;
  B.whr[2].re = 3;
  B.whr[3].re = 5;
  B.whr[4].re = 7;
  B.whr[5].re = 9;
  B.whr[6].re = 1;
  B.whr[7].re = 2;
  B.whr[8].re = 5;
    
  C.whr[0].re = -(14./5.+1./30.);
  C.whr[1].re = 2./3.;
  C.whr[2].re = .5;
  C.whr[3].re = 2+2./3.;
  C.whr[4].re = -1./3.;
  C.whr[5].re = -1;
  C.whr[6].re = -.5;
  C.whr[7].re = 0;
  C.whr[8].re = .5;

  B.whr[0].im = 0;
  B.whr[1].im = 0;
  B.whr[2].im = 0;
  B.whr[3].im = 0;
  B.whr[4].im = 0;
  B.whr[5].im = 0;
  B.whr[6].im = 0;
  B.whr[7].im = 0;
  B.whr[8].im = 0;
  
  C.whr[0].im = 0;
  C.whr[1].im = 0;
  C.whr[2].im = 0;
  C.whr[3].im = 0;
  C.whr[4].im = 0;
  C.whr[5].im = 0;
  C.whr[6].im = 0;
  C.whr[7].im = 0;
  C.whr[8].im = 0;

  fout << "\n\n#index = 2 - prodotto matrici scalari #\n" << endl;
  oo = (2 * 9 * 22);

  /* Calcolo fps: prodotto matrici scalari */
  for(int j = 0; j <= STOP-START; j++){
    data[j].flops  = 0;
    data[j].flops2 = 0;
    data[j].tempo  = 0;
    data[j].tempo2 = 0;

    for(int r = 0; r < RIPETIZIONI; r++){
      
      gettimeofday(&tv1, NULL);
      t1 = clock();
      for(int i = 0; i < data[j].i ; i++){
	A = A*B;
	A = A*C;
      }
      gettimeofday(&tv2, NULL);
      t2 = clock() -t1;
      printf("%f\n",A.whr[0].re);
      data[j].tempo  += t2/( double)CLOCKS_PER_SEC;
      data[j].tempo2 += ((tv2.tv_sec - tv1.tv_sec) +
			(tv2.tv_usec - tv1.tv_usec)/( double)1000000);
    }
    data[j].flops   = oo*RIPETIZIONI*data[j].i / data[j].tempo ;
    data[j].flops2  = oo*RIPETIZIONI*data[j].i / data[j].tempo2;
    
    fout << START+j         << "\t"
	 << scientific << data[j].tempo   << "\t"
	 << scientific << data[j].tempo2  << "\t"
	 << scientific << data[j].flops   << "\t"
	 << scientific << data[j].flops2  << endl;
  }


  /* Alloca matrici perturbative */
  ptSU3 CC;
  for(int i = 0; i < PTORD; i++){
    AA.ptU[i] = A;
    BB.ptU[i] = -A;
    CC.ptU[i] = A;
  }
  AA = exp(AA);
  BB = exp(BB);

  fout << "\n\n#index = 3 - prodotto matrici perturbative #\n" << endl;

  /* Calcolo fps: prodotto matrici perturbative */
  for(int j = 0; j <= STOP-START; j++){
    
    data[j].flops  = 0;
    data[j].flops2 = 0;
    data[j].tempo  = 0;
    data[j].tempo2 = 0;

    oo = ( 2 * (9 * PTORD*(PTORD+1) + 6 * (1+18*PTORD) + 11*9 *
		PTORD*(PTORD-1)) );
    for(int r = 0; r < RIPETIZIONI; r++){
      
      gettimeofday(&tv1, NULL);
      t1 = clock();
      for(int i = 0; i < data[j].i ; i++){
	CC = CC*AA;
	CC = CC*BB;
      }
      gettimeofday(&tv2, NULL);
      t2 = clock() - t1;
      printf("%f\n",CC.ptU[2].whr[0].re);
      data[j].tempo  += t2/( double)CLOCKS_PER_SEC;
      data[j].tempo2 += ((tv2.tv_sec - tv1.tv_sec) +
			 (tv2.tv_usec - tv1.tv_usec)/( double)1000000);
      
    }
    data[j].flops  = oo *data[j].i*RIPETIZIONI / data[j].tempo;
    data[j].flops2 = oo *data[j].i*RIPETIZIONI / data[j].tempo2;

    fout << START+j         << "\t"
	 << scientific << data[j].tempo   << "\t"
	 << scientific << data[j].tempo2  << "\t"
	 << scientific << data[j].flops   << "\t"
	 << scientific << data[j].flops2  << endl;
  }

  fout.close();

  system("gnuplot verifica.gp");
  return 0;
}
