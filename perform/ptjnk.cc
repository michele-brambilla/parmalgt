#include<iostream>
#include<fstream>
#include<stdio.h>
#include<QCDenvNODEpt.h>
#include<sys/time.h>


#define SUM_PT (2*9*PTORD)
#define SUM (2*9)
#define MULT (9*22)
#define MULT_C (6*9)
#define MULT_PT (.5*PTORD*(PTORD-1)*MULT+2*(PTORD-1)*MULT_C+PTORD*(PTORD+1)*SUM+6)
#define EXP_PT (PTORD*(PTORD*PTORD-1)/6*(MULT+SUM)+(PTORD-1)*(SUM_PT+MULT_C))

#define TG_INIT 2
#define TG_MAX 6

double ptjnk(int);



using namespace std;

int main(){

  long int seed = 124567;
#ifdef RANF
  rand_init(seed);
#else
  srand(seed);
#endif

  ofstream fout;
  fout.open("./result/ptjnk.dat");

  for(int i = TG_INIT; i <= TG_MAX; i += 2){
    fout << i << "\t" << ptjnk(i) << endl;
  }

  fout.close();
  return 0;

}

double ptjnk(int size){

  FILE *fin, *damocle;
  int SWEEP = 500;
  double tau = 0.02 , alpha = -.05;
  int sz[dim];

  unsigned long long operazioni;
  double flops;
  struct timeval tv1, tv2;
  clock_t t1, t2;
  long double tempo;
  
  sz[0] = sz[1] = sz[2] = sz[3] = size;
  
  latt LL(sz);  
  ptGluon_fld Umu(&LL);

  ptSU3 *U;
  int iVol = sz[0]*sz[1]*sz[2]*sz[3];

  ptSU3 F, FF;
  ptSU3 MM[dim], A, W1, W2, W;

  double stau = sqrt(tau), rVol;
  rVol = 1./(double)iVol;
  Cplx *w1, *ww1;
  w1  = new Cplx[PTORD+1];
  ww1 = new Cplx[PTORD+1];
  
  U = Umu.handle();

  for(int i = 0; i < iVol; i++){
    for(int mu = 0; mu < dim; mu++){
      U[get(&Umu, i, mu)].id();
    }
  }

  gettimeofday(&tv1, NULL);
  for(int t = 0; t < SWEEP; t++){
    
    for (int i1=0; i1 < dim; i1++){
      MM[i1].zero();
    }
    for (int i1=0; i1 <= PTORD; i1++){
      w1[i1] = 0.0;
    }
    for(int i = 0; i < iVol; i++){
      for(int mu = 0; mu < dim; mu++){
	
   	F.zero();
   	for(int nu = 0; nu < dim; nu++){
   	  if(nu != mu){
   	    F += Umu.staple(i, mu, nu);
	    // 3 staples = 
	    // 3*(4 m_p + 2 s_p) =
	  }
   	}
   	F = U[get(&Umu, i, mu)]*F;
	// 1 m_p
   	F.Tr(ww1);
   	for (int i1=0; i1 <= PTORD; i1++){
   	  w1[i1] += ww1[i1];
   	}
   	F *= -(tau/6.);
	// 1m_c
	F.reH();

   	F.ptU[0] -= stau*SU3rand();
	// 1m_c + 1s
	
   	U[get(&Umu, i, mu)] = exp(F)*U[get(&Umu, i, mu)];
	// 1exp + 1 m_p

	MM[mu] += log(U[get(&Umu, i, mu)]);
	// 1log + 1s_p
      }

    }
    // tot:
    // vol*4*(14 m_p + 7 s_p + 2 m_c + s + 2 ex +10o+44)
    // vol*(56 m_pt + 28 s_pt + 8 exp_pt + 8 m_c +4 s +40o +176)


//     // Zero momentum subtraction
//     for (int i1=0; i1 < dim; i1++){
//       MM[i1] *= rVol;
//     }
//     for(int i = 0; i < iVol; i++){ 
//       for(int mu = 0; mu < dim; mu++){
// 	A = log(U[get(&Umu, i, mu)]);
//  	A -= MM[mu];
// 	U[get(&Umu, i, mu)]  = exp(A);
//       }
//     }
//     //vol*( 8exp + 4s_p) + 4 m_c
    
//     // Gauge fixing
//     for(int i = 0; i < iVol; i++){
//       W.zero();
//       for(int mu = 0; mu < dim; mu++){
//  	W += ( U[get(&Umu, i, mu)] - dag(U[get(&Umu, i, mu)]) - 
//  	       U[get(&Umu, i, -1, mu, mu)] + 
// 	       dag(U[get(&Umu, i, -1, mu, mu)]) );
	
//  	W.Trless();
//       }
//       W1 = exp( alpha*W);
//       W2 = exp(-alpha*W);
//       for(int mu = 0; mu < dim; mu++){
//  	U[get(&Umu, i, mu)] = W1*U[get(&Umu, i, mu)];
//  	U[get(&Umu, i, -1, mu, mu)] = U[get(&Umu, i, -1, mu, mu)]*W2;
//       }
//     }
    // vol*( 16 s_pt + 16m_pt+8 exp_pt + 6o)

//     for (int i1=0; i1 <= PTORD; i1++){
//       fprintf(ff,"%le %le\n",w1[i1].re,w1[i1].im);
//     }    
    
    // load damocle.dag
//     if(t%BEAT == 0){
//       int sve, kll;
//       damocle = fopen("damocle.dag","r");
      
//       fscanf(damocle, "SWEEP %d \n",       &SWEEP);
//       fscanf(damocle, "BEAT %d \n",        &BEAT);
//       fscanf(damocle, "tau %lf\n",         &tau);
//       fscanf(damocle, "alpha %lf\n",       &alpha);
//       fscanf(damocle, "save_config %d \n", &sve);
//       fscanf(damocle, "kill_run %d   ",    &kll);    
      
//       fclose(damocle);
//       if(sve == 1) Umu.save(lastconf);
//       if(kll != 0) break;
//     }

  }
  gettimeofday(&tv2, NULL);
  tempo = (tv2.tv_sec - tv1.tv_sec) + 1e-6*(tv2.tv_usec - tv1.tv_usec);
  //  tempo = ( (tv2.tv_sec - tv1.tv_sec) + (tv2.tv_usec - tv1.tv_usec)/1000000.);

      

  operazioni = iVol*(56*MULT_PT + 28*SUM_PT + 8*EXP_PT + 8*MULT_C +4*SUM
		     +40*PTORD +176);
  

  flops = operazioni/tempo;
  cout << "Eseguite " << operazioni << "*SWEEP double operazioni in "
       << scientific
       << tempo << " secondi, per un effettivo di "
       << SWEEP*flops*1e-9 << " Gflops\n";

  delete []  w1;
  delete [] ww1;

  return (SWEEP*flops);
}
