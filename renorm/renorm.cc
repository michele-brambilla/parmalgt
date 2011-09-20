#include<RenLib.h>
#include<renorm_pars.h>
//#include<GammaStuff.h>

#include<fstream>
#include <sys/syscall.h>
#include <sys/types.h>
#include <unistd.h>
#include <omp.h>
#include <ctime>

int *sz, ptord, R, nthr, chunk;
char *conf, *out, *prfile, *prfileT;
double alpha;


using namespace std;

int main(int argc, char **argv){

  cout.precision(10);

  int cnt;

  sz = new int(dim);
  conf    = new char[100];
  out     = new char[100];
  prfile  = new char[100];
  prfileT = new char[100];

  if(get_params(argc, argv)){
    cout << "File di configurazione non trovato" << endl;
    //    exit(1);
  }

  latt LL(sz);
  LL.p_init();

  ptGluon_fld Umu(&LL);
  
  Cplx *plaq= new Cplx[allocORD+1];
  
  //  Umu.load_ape(conf);
  carica_campo(Umu,plaq);
  

  
  ptSU3 *U;
  U = Umu.handle();

  ptSU3 F,FF;
  int stat;

  if(gauge_fixing(Umu,cnt)){
    cout << "Attenzione, il gauge fixing non converge" << endl;
    exit(-1);
  }
  
  
  ptSpinColor_fld *Pmu[2*dim];
  
  ptSpinColor_fld Pmu0(&LL);
  ptSpinColor_fld Pmu1(&LL);
  ptSpinColor_fld Pmu2(&LL);
  ptSpinColor_fld Pmu3(&LL);
  ptSpinColor_fld Pmu4(&LL);
  ptSpinColor_fld Pmu5(&LL);
  ptSpinColor_fld Pmu6(&LL);
  ptSpinColor_fld Pmu7(&LL);

  Pmu[0] = &Pmu0;
  Pmu[1] = &Pmu1;
  Pmu[2] = &Pmu2;
  Pmu[3] = &Pmu3;
  Pmu[4] = &Pmu4;
  Pmu[5] = &Pmu5;
  Pmu[6] = &Pmu6;
  Pmu[7] = &Pmu7;


  Cplx propag[dim][dim][allocORD+1];

#ifdef _TENSOR_ 
  Cplx curr[16][dim][dim][NC][allocORD+1];
#else
  Cplx curr[10][dim][dim][NC][allocORD+1];
#endif



#ifdef _SAVE_TRANS_
  Cplx propagT[dim][dim][allocORD+1];
#endif

  ofstream momout;
  momout.open("running.log");

  FILE *fout, *pout;
#ifdef _SAVE_TRANS_
  FILE *tout;
#endif
  fout = fopen(out,"wb");
  pout = fopen(prfile,"wb");
#ifdef _SAVE_TRANS_
  tout = fopen(prfileT,"wb");
#endif

  ptSpinColor  gmp0[4],  gmp1[4],  gmp2[4],  gmp3[4], gmp4[4], gmp5[4];
  ptSpinColor gmp51[4], gmp52[4], gmp53[4], gmp54[4];
#ifdef _TENSOR_
  ptSpinColor gmp12[4], gmp13[4], gmp14[4];
  ptSpinColor gmp23[4], gmp24[4];
  ptSpinColor gmp34[4];
#endif
  ptSU3 scambio;


  int c;


  momout << "Configurazione: " << conf << endl;
  momout << "File Output Corrente: " << out << endl;
  momout << "File Propagatore: " << prfile << endl;
#ifdef _SAVE_TRANS_
  momout << "File PropagatoreTrasposto: " << prfileT << endl << endl;
#endif

  time_t tempo;
  time ( &tempo );
  momout << "Inizio run: " << ctime(&tempo) << endl << endl;
  momout << "Momenti eseguiti" << endl;


  FILE *fmom;
  if( (fmom = fopen("momREF.txt","r") ) == NULL){
	cout << "Errore, Momenti di Riferimento Non Trovati!" << endl;
	exit(1);
  }  
  int *xx       = new int[4];

  while( ( stat=fscanf(fmom,"%d %d %d %d\n", &xx[0], &xx[1], &xx[2], &xx[3])) >0){ ;
    
    c = LL.get(xx);
    cout << "x = \t"    <<  xx[0] << "\t"
	 << "y = \t"    <<  xx[1] << "\t"
	 << "z = \t"    <<  xx[2] << "\t"
	 << "t = \t"    <<  xx[3] << "\t"
	 << "sito = \t" <<  c     << "\t"
	 << "stato = "  <<  stat  << endl;
    
    	    
    memset(propag,0 ,sizeof(propag));
#ifdef _SAVE_TRANS_
    memset(propagT,0 ,sizeof(propagT));
#endif
    memset(curr,0 ,sizeof(curr));
    
    for(int a = 0; a < NC; a++){
      
      XiBuild1(Umu, (void**)Pmu, c, a);

#ifdef _OPEN_MP
#pragma omp barrier
#endif	      
      
      for(int mu = 0; mu < dim; mu++) {
	
	
	for(int o = 0 ; o <= ptord; o++){
	  
#ifdef _SAVE_TRANS_
	  propagT[0][mu][o] += Pmu[0]->psi[c].psi[mu].ptCV[o].whr[a];
	  propagT[1][mu][o] += Pmu[1]->psi[c].psi[mu].ptCV[o].whr[a];
	  propagT[2][mu][o] += Pmu[2]->psi[c].psi[mu].ptCV[o].whr[a];
	  propagT[3][mu][o] += Pmu[3]->psi[c].psi[mu].ptCV[o].whr[a];
#endif
	  
	  propag[mu][0][o] += Pmu[4]->psi[c].psi[mu].ptCV[o].whr[a];
	  propag[mu][1][o] += Pmu[5]->psi[c].psi[mu].ptCV[o].whr[a];
	  propag[mu][2][o] += Pmu[6]->psi[c].psi[mu].ptCV[o].whr[a];
	  propag[mu][3][o] += Pmu[7]->psi[c].psi[mu].ptCV[o].whr[a];
	  
	}
      }	      
      
      for( int d = 0; d < LL.Size; d++ ) {
	
	for(int i = 0; i < dim; i++)  {
	  gmp0[i]  = (Pmu[i+dim]->psi[d]);
	  gmp1[i]  = (Pmu[i+dim]->psi[d]).gmleft(0);
	  gmp2[i]  = (Pmu[i+dim]->psi[d]).gmleft(1);
	  gmp3[i]  = (Pmu[i+dim]->psi[d]).gmleft(2);
	  gmp4[i]  = (Pmu[i+dim]->psi[d]).gmleft(3);
	  gmp5[i]  = (Pmu[i+dim]->psi[d]).gmleft(5);
	  gmp51[i] = (gmp1[i]).gmleft(5);
	  gmp52[i] = (gmp2[i]).gmleft(5);
	  gmp53[i] = (gmp3[i]).gmleft(5);
	  gmp54[i] = (gmp4[i]).gmleft(5);
#ifdef _TENSOR_
	  gmp12[i] = (gmp1[i]).gmleft(5);
	  gmp13[i] = (gmp2[i]).gmleft(5);
	  gmp14[i] = (gmp3[i]).gmleft(5);
	  gmp23[i] = (gmp4[i]).gmleft(5);
	  gmp24[i] = (gmp4[i]).gmleft(5);
	  gmp34[i] = (gmp4[i]).gmleft(5);
#endif
	}
	
	
	for(int mu = 0; mu < dim; mu++) {
	  for(int nu = 0; nu < dim; nu++) {
	    
	    for( int iO = 0; iO <=ptord; iO++ ){
	      for( int jO = 0; jO <= iO; jO++ ){
		
		for( int rho = 0; rho < dim; rho++) {
		  for ( int icol = 0; icol < NC; icol++ ) {
		    
// #ifndef _CORRENTINE_
		    curr[0][mu][nu][a][iO] += 
		      Pmu[mu]->psi[d].psi[rho].ptCV[jO].whr[icol] * gmp0[nu].psi[rho].ptCV[iO-jO].whr[icol] ;
		    
		    curr[1][mu][nu][a][iO] += 
		      Pmu[mu]->psi[d].psi[rho].ptCV[jO].whr[icol] * gmp1[nu].psi[rho].ptCV[iO-jO].whr[icol] ; 
		    
		    curr[2][mu][nu][a][iO] += 
		      Pmu[mu]->psi[d].psi[rho].ptCV[jO].whr[icol] * gmp2[nu].psi[rho].ptCV[iO-jO].whr[icol] ; 

		    curr[3][mu][nu][a][iO] += 
		      Pmu[mu]->psi[d].psi[rho].ptCV[jO].whr[icol] * gmp3[nu].psi[rho].ptCV[iO-jO].whr[icol] ; 

		    curr[4][mu][nu][a][iO] += 
		      Pmu[mu]->psi[d].psi[rho].ptCV[jO].whr[icol] * gmp4[nu].psi[rho].ptCV[iO-jO].whr[icol] ; 

		    curr[5][mu][nu][a][iO] += 
		      Pmu[mu]->psi[d].psi[rho].ptCV[jO].whr[icol] * gmp5[nu].psi[rho].ptCV[iO-jO].whr[icol] ;

		    curr[6][mu][nu][a][iO] += 
		      Pmu[mu]->psi[d].psi[rho].ptCV[jO].whr[icol] * gmp51[nu].psi[rho].ptCV[iO-jO].whr[icol] ;

		    curr[7][mu][nu][a][iO] += 
		      Pmu[mu]->psi[d].psi[rho].ptCV[jO].whr[icol] * gmp52[nu].psi[rho].ptCV[iO-jO].whr[icol] ;

		    curr[8][mu][nu][a][iO] += 
		      Pmu[mu]->psi[d].psi[rho].ptCV[jO].whr[icol] * gmp53[nu].psi[rho].ptCV[iO-jO].whr[icol] ;

		    curr[9][mu][nu][a][iO] += 
		      Pmu[mu]->psi[d].psi[rho].ptCV[jO].whr[icol] * gmp54[nu].psi[rho].ptCV[iO-jO].whr[icol] ;
// #else
// 		    curr[0][mu][nu][a][iO] += 
// 		      Pmu[mu]->psi[d].psi[rho].ptCV[jO].whr[icol] * gmp51[nu].psi[rho].ptCV[iO-jO].whr[icol] ;

// 		    curr[1][mu][nu][a][iO] += 
// 		      Pmu[mu]->psi[d].psi[rho].ptCV[jO].whr[icol] * gmp52[nu].psi[rho].ptCV[iO-jO].whr[icol] ;

// 		    curr[2][mu][nu][a][iO] += 
// 		      Pmu[mu]->psi[d].psi[rho].ptCV[jO].whr[icol] * gmp53[nu].psi[rho].ptCV[iO-jO].whr[icol] ;

// 		    curr[3][mu][nu][a][iO] += 
// 		      Pmu[mu]->psi[d].psi[rho].ptCV[jO].whr[icol] * gmp54[nu].psi[rho].ptCV[iO-jO].whr[icol] ;
// #endif
		  }
		}
	      }
	    }
	  }
	}		
      }
	      
    }
	        
    fwrite(curr,    sizeof(curr),    1, fout);
    fwrite(propag,  sizeof(propag),  1, pout);
#ifdef _SAVE_TRANS_
    fwrite(propagT, sizeof(propagT), 1, tout);
#endif

    momout << xx[0] << "\t"
	   << xx[1] << "\t"
	   << xx[2] << "\t"
	   << xx[3] << "\n";
    momout.flush();
    
  }
  fclose(fmom);
  
  time ( &tempo );
  momout << "Fine run: " << ctime(&tempo) << endl << endl;
    
  
  //   delete Pmu;
  fclose(fout);
  fclose(pout);
#ifdef _SAVE_TRANS_
  fclose(tout);
#endif
  momout.close();
  return 0;
}

