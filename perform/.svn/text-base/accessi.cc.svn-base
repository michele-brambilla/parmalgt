/* Statistica sui tempi di accesso utilizzando la classe reticolo */

#include<iostream>
#include<fstream>
#include<stdio.h>
#include<QCDenvNODEpt.h>
#include<sys/time.h>

#define SWEEP 100000  // numero massimo cicli
#define RIP 5         // ripetizioni per statistica
#define TG_INIT 2      // taglia minima
#define TG_MAX 8       // taglia massima

using namespace std;

// riceve taglia, file output e tempo di copia
void accessi_metodo(int, ofstream&, double);
void accessi_friend(int, ofstream&, double);


int main(){
  long int seed = 1234567;
#ifdef RANF
  rand_init(seed);
#else
  srand(seed);
#endif

  SU3 A,B;
  A = SU3rand();
  cout << "Calibro velocita' copia..." << endl;
  timeval t1,t2, t;
  double tmp;
  t.tv_sec  = 0;
  t.tv_usec = 0;

  for(int j = 0; j < C_RP; j++){
    gettimeofday(&t1, NULL);
    for(int i = 0; i < C_SW; i++){
      B = A;
    }
    gettimeofday(&t2, NULL);
    t.tv_sec  += t2.tv_sec  - t1.tv_sec; 
    t.tv_usec += t2.tv_usec - t1.tv_usec;
  }
  cout << (tmp = (t.tv_sec*1000000. + 
		 t.tv_usec)/(double)(C_SW*C_RP) )
       << "\tmsec" << endl;
  cout << "Pronto" << endl;

  int sz;
  cout.precision(5);

  ofstream fout;
  fout.open("./result/accessi_m.dat");
  fout << "#it\tgettimeofday()" << endl;
  for(int tg = TG_INIT; tg <= TG_MAX; tg += 2){
    fout << "# taglia = " << tg << endl;
    accessi_metodo(tg, fout, tmp*1e-6);
    fout << "\n\n";
  }
  fout.close();

  fout.open("./result/accessi_f.dat");
  fout << "#it\tgettimeofday()" << endl;
  for(int tg = TG_INIT; tg <= TG_MAX; tg += 2){
    fout << "# taglia = " << tg << endl;
    accessi_friend(tg, fout, tmp*1e-6);
    fout << "\n\n";
  }
  fout.close();

  system("gnuplot accessi.gp");
  return 0;
}



void accessi_friend(int size, ofstream &fout, double tmp){
  struct timeval tv1, tv2;
  clock_t tempo;

  int sz[4];
  sz[0] = sz[1] = sz[2] = sz[3] = size;

  latt LL(sz);  
  ptGluon_fld Umu(&LL);

  ptSU3 *U; 
  SU3 A, B;
  int iVol = sz[0]*sz[1]*sz[2]*sz[3] , norm;

  U = Umu.handle();

  for(int i = 0; i < iVol; i++){
    for(int mu = 0; mu < dim; mu++){
      U[get(&Umu, i, mu)].id();
    }
  }
      
  for(long int max = 10; max < SWEEP; max*=10){
    norm = max*RIP*pow((double)size,4);
    tempo = 0;
    for(int r = 0; r < RIP; r++){
      gettimeofday(&tv1, NULL);
      for(long int t = 0; t < max; t++){
	
	for(int i = 0; i < iVol; i++){
	  for(int mu = 0; mu < dim; mu++){
	    A = U[get(&Umu, i, mu)].ptU[0];
	  }
	}
      }

      gettimeofday(&tv2, NULL);
      tempo += ( (tv2.tv_sec - tv1.tv_sec)*1000000. + 
		 (tv2.tv_usec - tv1.tv_usec));
      A.prout();
    }
    fout << (int)log10(max)
	 << "\t"
	 << scientific << tempo*1e-6/(double)(norm) - tmp << "\n";
    fout.flush();
    cout << (int)log10(max) << endl;
    fflush(stdout);
  }
}


void accessi_metodo(int size, ofstream &fout, double tmp){
  struct timeval tv1, tv2;
  clock_t tempo;
  int sz[4];
  sz[0] = sz[1] = sz[2] = sz[3] = size;

  latt LL(sz);  
  ptGluon_fld Umu(&LL);

  ptSU3 *U; 
  SU3 A, B;
  int iVol = sz[0]*sz[1]*sz[2]*sz[3] , norm;

  U = Umu.handle();

  for(int i = 0; i < iVol; i++){
    for(int mu = 0; mu < dim; mu++){
      U[get(&Umu, i, mu)].id();
    }
  }
      
  for(long int max = 10; max < SWEEP; max*=10){
    norm = max*RIP*pow((double)size,4);
    tempo = 0;
    for(int r = 0; r < RIP; r++){
      gettimeofday(&tv1, NULL);
      for(long int t = 0; t < max; t++){
	for(int i = 0; i < iVol; i++){
	  for(int mu = 0; mu < dim; mu++){
	    A = (Umu.get(i, mu)).ptU[0];
	  }
	}
      }

      gettimeofday(&tv2, NULL);
      tempo += ( (tv2.tv_sec - tv1.tv_sec)*1000000. + 
		 (tv2.tv_usec - tv1.tv_usec));
      A.prout();
    }
    fout << (int)log10(max)
	 << "\t"
	 << scientific << tempo*1e-6/(double)(norm) - tmp << "\n";
    fout.flush();
    cout << (int)log10(max) << endl;
    fflush(stdout);
  }
}
