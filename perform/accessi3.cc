#include<iostream>
#include<fstream>
#include<stdio.h>
#include<QCDenvNODEpt.h>
#include<sys/time.h>

#define SWEEP 10       // numero massimo cicli
#define RIP 10         // ripetizioni per statistica
#define TG_INIT 4      // taglia minima
#define TG_MAX 8       // taglia massima

using namespace std;

// riceve taglia, file output e tempo di copia
timeval accessi_metodo(int);
timeval accessi_friend(int);


int main(){
  long int seed = 1234567;
#ifdef RANF
  rand_init(seed);
#else
  srand(seed);
#endif

  SU3 A,B;
  timeval t1,t2, tf, tm, t;
  //  timev t;
  long double norm;
  ofstream fout, mout;
  mout.open("./result/accessi_m3.dat");
  fout.open("./result/accessi_f3.dat");

  cout.precision(10);
  
  mout << "# taglia\tcalibrazione\toperazione\taccesso"<< endl;
  fout << "# taglia\tcalibrazione\toperazione\taccesso"<< endl;
  

  for(int tg = TG_INIT; tg <= TG_MAX; tg += 2){

    A = 0.00000001*SU3rand();
    A.whr[0] +=1;
    A.whr[4] +=1;
    A.whr[8] +=1;
    B = A;
    B.prout();
    t.tv_sec  = 0;
    t.tv_usec = 0;
    int iVolMax = tg*tg*tg*tg;

    cout << "Calibro velocita' operazioni: taglia " << tg << endl;
    
    for(int j = 0; j < SWEEP; j++){
      gettimeofday(&t1, NULL);
      for(int i = 0; i < RIP; i++){
	for(int v = 0; v < iVolMax; v++){
	  for(int mu = 0; mu < dim; mu++){
	    B *= A;
	    B *= A;
	    B *= A;
	    B *= A;
	    B *= A;
	  }
	}
      }
      gettimeofday(&t2, NULL);
      t.tv_sec  += t2.tv_sec  - t1.tv_sec; 
      t.tv_usec += t2.tv_usec - t1.tv_usec;
//       cout << scientific
// 	   << t2.tv_sec << "\t" << t1.tv_sec << "\t" 
// 	   << t2.tv_usec   << "\t" <<  t1.tv_usec << "\t"
// 	   << (t2.tv_sec - t1.tv_sec) + 1e-6*(t2.tv_usec - t1.tv_usec)
// 	   << endl;
    }
    cout << scientific << "\n"
	 << (t.tv_sec + 1e-6*t.tv_usec)/(long double)(SWEEP*dim*5*RIP*iVolMax) 
	 << endl;

    cout << "Pronto" << endl;
    B.prout();

    norm = 1./(iVolMax*SWEEP*RIP*5);

    tm = accessi_metodo(tg);
    tf = accessi_friend(tg);
    
    mout << tg << "\t"
	 << scientific 
	 << t.tv_sec  + 1e-6*t.tv_usec  << "\t"
	 << tm.tv_sec + 1e-6*tm.tv_usec << "\t"
	 << ((tm.tv_sec-t.tv_sec) + 1e-6*(tm.tv_usec-t.tv_usec))*norm
	 << endl;

    fout << tg << "\t"
	 << scientific 
	 << t.tv_sec  + 1e-6*t.tv_usec  << "\t"
	 << tf.tv_sec + 1e-6*tf.tv_usec << "\t"
	 << ((tf.tv_sec-t.tv_sec)+1e-6*(tf.tv_usec-t.tv_usec))*norm
	 << endl;

    mout.flush();
    fout.flush();
  }



  fout.close();
  mout.close();

  return 0;
}



timeval accessi_friend(int size){
  struct timeval tv1, tv2, tempo;

  int sz[4];
  sz[0] = sz[1] = sz[2] = sz[3] = size;

  latt LL(sz);  
  ptGluon_fld Umu(&LL);

  ptSU3 *U; 
  SU3 A, B, C, D, E, F;
  int iVol = sz[0]*sz[1]*sz[2]*sz[3] , norm;

  A = 0.00000001*SU3rand();
  A.whr[0] +=1;
  A.whr[4] +=1;
  A.whr[8] +=1;
  B = A;
  C = A;

  U = Umu.handle();

  for(int i = 0; i < iVol; i++){
    for(int mu = 0; mu < dim; mu++){
      U[get(&Umu, i, mu)].id();
      U[get(&Umu, i, mu)].ptU[0] = 0.00000001*SU3rand();
      U[get(&Umu, i, mu)].ptU[0].whr[0] +=1;
      U[get(&Umu, i, mu)].ptU[0].whr[4] +=1;
      U[get(&Umu, i, mu)].ptU[0].whr[8] +=1;
      U[get(&Umu, i, mu)].ptU[1] = 0.00000001*SU3rand();
      U[get(&Umu, i, mu)].ptU[1].whr[0] +=1;
      U[get(&Umu, i, mu)].ptU[1].whr[4] +=1;
      U[get(&Umu, i, mu)].ptU[1].whr[8] +=1;
      U[get(&Umu, i, mu)].ptU[2] = 0.00000001*SU3rand();
      U[get(&Umu, i, mu)].ptU[2].whr[0] +=1;
      U[get(&Umu, i, mu)].ptU[2].whr[4] +=1;
      U[get(&Umu, i, mu)].ptU[2].whr[8] +=1;
      U[get(&Umu, i, mu)].ptU[3] = 0.00000001*SU3rand();
      U[get(&Umu, i, mu)].ptU[3].whr[0] +=1;
      U[get(&Umu, i, mu)].ptU[3].whr[4] +=1;
      U[get(&Umu, i, mu)].ptU[3].whr[8] +=1;
      U[get(&Umu, i, mu)].ptU[4] = 0.00000001*SU3rand();
      U[get(&Umu, i, mu)].ptU[4].whr[0] +=1;
      U[get(&Umu, i, mu)].ptU[4].whr[4] +=1;
      U[get(&Umu, i, mu)].ptU[4].whr[8] +=1;
    }
  }
      
  tempo.tv_sec = tempo.tv_usec = 0;

  for(long int t = 0; t < SWEEP; t++){
    gettimeofday(&tv1, NULL);
    for(int r = 0; r < RIP; r++){
      
      for(int i = 0; i < iVol; i++){
	for(int mu = 0; mu < dim; mu++){
	  C = U[get(&Umu, i, mu)].ptU[0];
	  A *= B;
	  D = U[get(&Umu, i, mu)].ptU[1];
	  A *= C;
	  E = U[get(&Umu, i, mu)].ptU[2];
	  A *= D;
	  F = U[get(&Umu, i, mu)].ptU[3];
	  A *= E;
	  B = U[get(&Umu, i, mu)].ptU[4];
	  A *= F;
	}
      }
    }
    
    gettimeofday(&tv2, NULL);

    tempo.tv_sec  += tv2.tv_sec - tv1.tv_sec;
    tempo.tv_usec += tv2.tv_usec - tv1.tv_usec;

    A.prout();
  }

  return tempo;
}


timeval accessi_metodo(int size){
  struct timeval tv1, tv2, tempo;
  int sz[4];
  sz[0] = sz[1] = sz[2] = sz[3] = size;

  latt LL(sz);  
  ptGluon_fld Umu(&LL);

  ptSU3 *U; 
  SU3 A, B, C, D, E, F;
  int iVol = sz[0]*sz[1]*sz[2]*sz[3] , norm;

  A = 0.00000001*SU3rand();
  A.whr[0] +=1;
  A.whr[4] +=1;
  A.whr[8] +=1;
  B = A;
  C = A;

  U = Umu.handle();

  for(int i = 0; i < iVol; i++){
    for(int mu = 0; mu < dim; mu++){
      U[get(&Umu, i, mu)].id();
      U[get(&Umu, i, mu)].ptU[0] = 0.00000001*SU3rand();
      U[get(&Umu, i, mu)].ptU[0].whr[0] +=1;
      U[get(&Umu, i, mu)].ptU[0].whr[4] +=1;
      U[get(&Umu, i, mu)].ptU[0].whr[8] +=1;
      U[get(&Umu, i, mu)].ptU[1] = 0.00000001*SU3rand();
      U[get(&Umu, i, mu)].ptU[1].whr[0] +=1;
      U[get(&Umu, i, mu)].ptU[1].whr[4] +=1;
      U[get(&Umu, i, mu)].ptU[1].whr[8] +=1;
      U[get(&Umu, i, mu)].ptU[2] = 0.00000001*SU3rand();
      U[get(&Umu, i, mu)].ptU[2].whr[0] +=1;
      U[get(&Umu, i, mu)].ptU[2].whr[4] +=1;
      U[get(&Umu, i, mu)].ptU[2].whr[8] +=1;
      U[get(&Umu, i, mu)].ptU[3] = 0.00000001*SU3rand();
      U[get(&Umu, i, mu)].ptU[3].whr[0] +=1;
      U[get(&Umu, i, mu)].ptU[3].whr[4] +=1;
      U[get(&Umu, i, mu)].ptU[3].whr[8] +=1;
      U[get(&Umu, i, mu)].ptU[4] = 0.00000001*SU3rand();
      U[get(&Umu, i, mu)].ptU[4].whr[0] +=1;
      U[get(&Umu, i, mu)].ptU[4].whr[4] +=1;
      U[get(&Umu, i, mu)].ptU[4].whr[8] +=1;
    }
  }
      
  tempo.tv_sec = tempo.tv_usec = 0;

  for(long int t = 0; t < SWEEP; t++){
    gettimeofday(&tv1, NULL);
    for(int r = 0; r < RIP; r++){
      for(int i = 0; i < iVol; i++){
	for(int mu = 0; mu < dim; mu++){
	  C = Umu.get(i, mu).ptU[0];
	  A *= B;
	  D = Umu.get(i, mu).ptU[1];
	  A *= C;
	  E = Umu.get(i, mu).ptU[2];
	  A *= D;
	  F = Umu.get(i, mu).ptU[3];
	  A *= E;
	  B = Umu.get(i, mu).ptU[4];
	  A *= F;
	}
      }
    }
    
    gettimeofday(&tv2, NULL);

    tempo.tv_sec  += tv2.tv_sec - tv1.tv_sec;
    tempo.tv_usec += tv2.tv_usec - tv1.tv_usec;

    A.prout();
  }

  return tempo;
}
