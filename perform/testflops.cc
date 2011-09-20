#include<iostream>
#include<fstream>
#include<stdio.h>
#include"QCDenvNODEpt.h"
#include<sys/time.h>

#define MULT_PT (9*PTORD*(PTORD+1)+6*(1+18*PTORD)+11*9*PTORD*(PTORD-1))
#define SUM_PT (6*9*PTORD)
#define SUM (6*9)
#define MULT (9*22)
#define MULT_C (6*9)
#define EXP_PT (PTORD+1)*PTORD*(PTORD-1)/6*MULT_PT + ((PTORD+1)*PTORD*(PTORD-1)/6 + PTORD-1)*SUM_PT



using namespace std;

int PTjnk(clock_t &tv, clock_t &tmp, int taglia){
  FILE *fin, *damocle;
  int SWEEP, BEAT;
  long int seed;
  char  str[100], lastconf[100];
  double tau, alpha;
  bool status;
  int sz[dim],term;
  unsigned long long operazioni, flops;
  struct timeval tv1, tv2;
  clock_t t1, t2, tempo;

  fin = fopen("ptconf.txt","r");

  fscanf(fin, "taglia %d %d %d %d \n", &sz[0], &sz[1], &sz[2], &sz[3]);
  fscanf(fin, "SWEEP %d \n",           &SWEEP);
  fscanf(fin, "BEAT %d \n",            &BEAT);
  fscanf(fin, "seed %ld\n",            &seed);
  fscanf(fin, "tau %lf\n",             &tau);
  fscanf(fin, "alpha %lf\n",           &alpha);
  fscanf(fin, "init_status %d \n",     &status);
  fscanf(fin, "no_term %d \n",         &term);
  fscanf(fin, "plaq_out %s \n",        str);
  fscanf(fin, "last_conf %s   ",       lastconf);
  fclose(fin);
  sz[0] = sz[1] = sz[2] = sz[3] = taglia;

#ifdef RANF
  rand_init(seed);
#else
  srand(seed);
#endif

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
  FILE *ff;
  
  U = Umu.handle();

  if(!status){
    for(int i = 0; i < iVol; i++){
      for(int mu = 0; mu < dim; mu++){
	U[get(&Umu, i, mu)].id();
      }
    }
  }
  else{
    Umu.load(lastconf);
  }
      
  gettimeofday(&tv1, NULL);
  t1 = clock();
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
    // vol*4*(3*(4 m_p + 2 s_p) + 2m_p + 2m_c + 1 s_p + 1s +2exp)
    // vol*4*(14 m_p + 7 s_p + 2 m_c + s + 2 ex)


    // Zero momentum subtraction
    for (int i1=0; i1 < dim; i1++){
      MM[i1] *= rVol;
    }
    for(int i = 0; i < iVol; i++){ 
      for(int mu = 0; mu < dim; mu++){
	A = log(U[get(&Umu, i, mu)]);
 	A -= MM[mu];
	U[get(&Umu, i, mu)]  = exp(A);
      }
    }
    //vol*4*( exp + s_p)
    
    // Gauge fixing
    for(int i = 0; i < iVol; i++){
      W.zero();
      for(int mu = 0; mu < dim; mu++){
 	W += ( U[get(&Umu, i, mu)] - dag(U[get(&Umu, i, mu)]) - 
 	       U[get(&Umu, i, -1, mu, mu)] + 
	       dag(U[get(&Umu, i, -1, mu, mu)]) );
	
 	W.Trless();
      }
      W1 = exp( alpha*W);
      W2 = exp(-alpha*W);
      for(int mu = 0; mu < dim; mu++){
 	U[get(&Umu, i, mu)] = W1*U[get(&Umu, i, mu)];
 	U[get(&Umu, i, -1, mu, mu)] = U[get(&Umu, i, -1, mu, mu)]*W2;
      }
    }
    // vol*( 2 exp + 2 m_p + 16s_p)

  }
  t2 = clock() - t1;
  gettimeofday(&tv2, NULL);

  tv = ( (tv2.tv_sec - tv1.tv_sec) + 
	  (tv2.tv_usec - tv1.tv_usec)*1e-6);
  tmp = t2;

  // tot:
  //  operazioni = 2*iVol*(9*EXP+77*MULT_PT+38*SUM_PT+4*MULT_C); 
  // 2 (216 + 9 ex + 2052 PTORD + 77 (99 (-1 + PTORD) PTORD + 
  // 9 PTORD (1 + PTORD) + 6 (1 + 18 PTORD))) vol
  // 4^4 PTORD = 8  --> 3083642880 fpo
  delete []  w1;
  delete [] ww1;
  return 0;
}

#define RIPETIZIONI 4
#define N_TAGLIE 3

int main(){
  clock_t tv, tempi, ts, tmps;
  int taglia[N_TAGLIE] = {4,6,8};
  ofstream outf;
  outf.open("./result/prlgt_time.dat");

  outf << "# taglia\tclock\tgettimeofday" << endl;
  for(int i = 0; i < N_TAGLIE; i++){
    ts = 0.;
    tmps = 0.;
    for(int rip = 0; rip < RIPETIZIONI; rip++){
      PTjnk(tv, tempi, taglia[i]);
      ts += tv;
      tmps += tempi;
    }
    cout << "Taglia: " << taglia[i]
	 << "\tClock: " << tmps/(double)(RIPETIZIONI*pow(taglia[i],4))
	 << "\tGettimeofday: " << ts/(double)(RIPETIZIONI*pow(taglia[i],4))
	 << endl;

    outf << taglia[i] << "\t"
	 << tmps/(double)RIPETIZIONI << "\t"
	 << ts/(double)RIPETIZIONI
	 << endl;
  }

  outf.close();
  return 0;
}
