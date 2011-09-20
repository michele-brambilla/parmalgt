#include "nspt.h"

ptSU2 W, W1, W2, Wsym;
ptSU2 *MM = new ptSU2[dim];

int link_c, site_c, site_up, xx[4], curr;
const int offset = dim+1;

extern Cplx* w;  // usero' come variabile d'appoggio
extern Cplx* w1; // usero' per placchetta 1x1
extern Cplx* w2; // usero' per placchetta 2x1
extern double* norm;
extern double* norm1;
extern Cplx** trU;
extern int ratio; 
extern ptSU2* U;
extern nspt_params_t nspt_pars;
extern act_params_t  act_pars;
extern thread_params_t thr_pars[NTHR];

extern std::ofstream plaqfile, normfile, logfile, trufile;

extern PRlgtTime Time;

#ifdef __WILLOOP_2x2__
extern std::ofstream loop2x2;
#endif

#ifndef __PARALLEL_OMP__
extern MyRand Rand;
#else
ptSU2 Ww[NTHR], Ww1[NTHR], Ww2[NTHR], Wwsym[NTHR];
ptSU2 MMm[NTHR][dim];
SU2 P[NTHR];
Cplx ww[NTHR][allocORD+1], ww1[NTHR][allocORD+1], ww2[NTHR][allocORD+1];
Cplx ttrU[NTHR][dim][allocORD];
double normm[NTHR][allocORD];
int tid;

extern MyRand Rand[];

#endif


void gauge_wilson(ptBoson_fld& Umu){

#ifndef __PARALLEL_OMP__

  memset(norm, 0, PTORD*sizeof(double));

  for(curr = 0; curr < act_pars.iVol; curr++){
    
    for(int mu = 0; mu < dim; mu++){
      // Preparo il link corrente
      link_c = get(&Umu, curr, mu);
      
      W1.zero();	
#if GAUGE_ACTION != WIL
      W2.zero();
#endif				
      // Calcolo la staple, eventualmente alnche la 2x1 su tutte
      // le direzioni ortogonali al link corrente
      for(int nu = 0; nu < dim; nu++){		
	if(nu != mu){					
	  W1 += Umu.staple(curr, mu, nu);			
#if GAUGE_ACTION != WIL
	  W2 += Umu.staple2x1(curr, mu, nu);			
#endif
	}						
      }						
      
      // Chido le staple a dare la placchetta e traccio
      W1 = U[link_c]*W1;
      W1.Tr(w);
      // Sommo al valore della placchetta e calcolo la norma
      w1[0] += w[0];
      for (int i1=1; i1 <= PTORD; i1++){		
	w1[i1] += w[i1];
	norm[i1-1] += ((U[link_c].ptU[i1-1]*dag(U[link_c].ptU[i1-1])).Tr()).re;
	trU[mu][i1-1] += U[link_c].ptU[i1-1].Tr();
      }
      
      // Se ho azioni improved faccio lo stesso sulle
      // parti di improvement
#if GAUGE_ACTION != WIL
      W2 = U[link_c]*W2;			
      W2.Tr(w);
      for (int i1=0; i1 <= PTORD; i1++){		
	w2[i1] += w[i1];
      }
      
      W1 = act_pars.c0*W1 + act_pars.c1*W2;
#endif
      // Costringo nell'algebra
      W1 = W1.reH();					
      W1 *= act_pars.tau_g;
      
      // Aggiungo fluttuazione al primo ordine
      W1.ptU[0] -= act_pars.stau*SU2rand(Rand);
      U[link_c] = exp(W1)*U[link_c];	
      
      // Contributo del momento nullo
      MM[mu] += log(U[link_c]);		
      
      
    } //end mu
    
  } // loop sui siti


#else  // if parallel

  for(tid = 0;tid < NTHR; tid++) {
    for(int ord = 0; ord < PTORD; ord++){
	ww[tid][ord]  = 0;
	ww1[tid][ord] = 0;
	ww2[tid][ord] = 0;
	normm[tid][ord] = 0;
	for( int mu = 0; mu < dim; mu++)
	  {
	    ttrU[tid][mu][ord] = Cplx(0,0);
	  }
    }
    ww[tid][allocORD]  = 0;
    ww1[tid][allocORD] = 0;
    ww2[tid][allocORD] = 0;
    
    
    for(int mu = 0; mu < dim; mu++){
      MMm[tid][mu].zero();
    }
  }
  


#pragma omp parallel private(curr,tid,link_c) num_threads(NTHR) 
  {
    tid = omp_get_thread_num();
    for(int y0 = thr_pars[tid].xi[0]; y0 < thr_pars[tid].xf[0]; y0++){
      for(int y1 = thr_pars[tid].xi[1]; y1 < thr_pars[tid].xf[1]; y1++){
	for(int y2 = thr_pars[tid].xi[2]; y2 < thr_pars[tid].xf[2]; y2++){
	  for(int y3 = thr_pars[tid].xi[3]; y3 < thr_pars[tid].xf[3]; y3++){

	    curr = y3 + act_pars.sz[3]*(y2 + act_pars.sz[2]*(y1 + act_pars.sz[1]*y0) );


	    for(int mu = 0; mu < dim; mu++){
	      // Preparo il link corrente
	      link_c = get(&Umu, curr, mu);
      
	      Ww1[tid].zero();	
#if GAUGE_ACTION != WIL
	      Ww2[tid].zero();
#endif				
	      // Calcolo la staple, eventualmente anche la 2x1 su tutte
	      // le direzioni ortogonali al link corrente
	      for(int nu = 0; nu < dim; nu++){		
		if(nu != mu){					
		  Ww1[tid] += Umu.staple(curr, mu, nu);			
#if GAUGE_ACTION != WIL
		  Ww2[tid] += Umu.staple2x1(curr, mu, nu);			
#endif
		}						
	      }						
      
	      // Chido le staple a dare la placchetta e traccio
	      Ww1[tid] = U[link_c]*Ww1[tid];
	      Ww1[tid].Tr(ww[tid]);
	      // Sommo al valore della placchetta e calcolo la norma
	      ww1[tid][0] += ww[tid][0];
	      for (int i1=1; i1 <= PTORD; i1++){
		//std::cout << "\t\t\t\ti1 =" << i1 << std::endl;
			      
		ww1[tid][i1] += ww[tid][i1];
		normm[tid][i1-1] += ((U[link_c].ptU[i1-1]*dag(U[link_c].ptU[i1-1])).Tr()).re;

		ttrU[tid][mu][i1-1] += U[link_c].ptU[i1-1].Tr();
	      }
      
	      // Se ho azioni improved faccio lo stesso sulle
	      // parti di improvement
#if GAUGE_ACTION != WIL
	      Ww2[tid] = U[link_c]*Ww2[tid];			
	      Ww2[tid].Tr(ww[tid]);
	      for (int i1=0; i1 <= PTORD; i1++){		
		ww2[tid][i1] += ww[tid][i1];
	      }
      
	      Ww1[tid] = act_pars.c0*Ww1[tid] + act_pars.c1*Ww2[tid];
#endif
	      // Costringo nell'algebra
	      Ww1[tid]  = Ww1[tid].reH();					
	      Ww1[tid] *= act_pars.tau_g;
      
	      // Aggiungo fluttuazione al primo ordine
	      P[tid] = SU2rand(Rand[tid]);
	      Ww1[tid].ptU[0] -= act_pars.stau*P[tid];
	      U[link_c] = exp(Ww1[tid])*U[link_c];	
      
	      // Contributo del momento nullo
	      MMm[tid][mu] += log(U[link_c]);		
      
	    } //end mu
    
#if ntz > 1
#pragma omp barrier
#endif
	  } // end z
	  
#if nty > 1
#pragma omp barrier
#endif
	} // end y
	
#if ntx > 1
#pragma omp barrier
#endif
      } // end x
      
#if ntt > 1
#pragma omp barrier
#endif
    }  // end t

#pragma omp barrier  
  } // end parallel


  // Reduce from parallel
  for(tid = 0; tid < NTHR; tid++){
  
    for(int ord = 0; ord < PTORD; ord++){
      w1[ord] += ww1[tid][ord];
      w2[ord] += ww2[tid][ord];
      for(int mu = 0; mu < dim; mu++){
	trU[mu][ord] += ttrU[tid][mu][ord];
      }
      norm[ord] += normm[tid][ord];
    }
    w1[PTORD] += ww1[tid][PTORD];
    w2[PTORD] += ww2[tid][PTORD];

    for(int mu = 0; mu < dim; mu++){
      MM[mu] += MMm[tid][mu];
    }
  }

#endif


  for (int i1=0; i1 < PTORD; i1++){		
#if GAUGE_ACTION == WIL
    plaqfile << w1[i1].re/(double)(48*act_pars.iVol) << "\t" 
	     << w1[i1].im/(double)(48*act_pars.iVol) << std::endl;
#else
    plaqfile << w1[i1].re/(double)(48*act_pars.iVol) << "\t" 
	     << w1[i1].im/(double)(48*act_pars.iVol) << "\t"
	     << w2[i1].re/(double)(48*act_pars.iVol) << "\t" 
	     << w2[i1].im/(double)(48*act_pars.iVol) << std::endl;
#endif
    normfile << norm[i1]/(double)(2*dim*act_pars.iVol) << std::endl;
    for( int mu = 0; mu < dim; mu++)
      {
	trufile  << trU[mu][i1].re /(double)(3*act_pars.iVol) << "\t" 
		 << trU[mu][i1].im /(double)(3*act_pars.iVol) << "\t";
      }
    trufile  << std::endl;
  }

#if GAUGE_ACTION == WIL
  plaqfile << w1[PTORD].re/(double)(48*act_pars.iVol) << "\t" 
	   << w1[PTORD].im/(double)(48*act_pars.iVol) << std::endl;
#else
  plaqfile << w1[PTORD].re/(double)(48*act_pars.iVol) << "\t" 
	   << w1[PTORD].im/(double)(48*act_pars.iVol) << "\t"
	   << w2[PTORD].re/(double)(48*act_pars.iVol) << "\t" 
	   << w2[PTORD].im/(double)(48*act_pars.iVol) << std::endl;
#endif

  
} // fine evoluzione gluoni




// SpinColor ZeroMom;
// SpinColor XiTmp;
// ptSpinColor PsiTmp;

// void fermion_wilson(ptBoson_fld& Umu, ptSpinColor_fld& Pmu, SpinColor_fld& Xi, ptBoson_fld& Wfld){

//   // genero campo gaussiano
//   Xi.gauss();

//   // Sottraggo il momento nullo
//   memset((void*)&ZeroMom, 0, sizeof(SpinColor));
//   for( int iv = 0; iv < act_pars.iVol; iv++){
//     ZeroMom += Xi.psi[iv];
//   }
//   ZeroMom *= act_pars.rVol;
//   for( int iv = 0; iv < act_pars.iVol; iv++){
//     Xi.psi[iv] -= ZeroMom;
//   }

//   // inverto sul campo generato l'operatore di Dirac.
//   // Prima copio il campo Xi in un campo ausiliario contenuto
//   // nel campo perturbativo
//   (*Pmu.scfld) = Xi;
 
//   memset((void*)Pmu.psi, 0, act_pars.iVol*sizeof(ptSpinColor));
//   Pmu.fillPT(PTORD-2,Umu);
  
//   // cicla su tutti i siti
//   for(curr = 0; curr < act_pars.iVol; curr++){
//     site_c = Umu.Z->L[curr][dim];

//     for( int mu = 0; mu < dim; mu++ ){
      
//       // variabile d'appoggio
//       W.zero();
//       XiTmp   = Xi.psi[Umu.Z->L[curr][offset+mu]];
//       PsiTmp  = Pmu.psi[site_c];
//       PsiTmp += PsiTmp.gmleft(mu);
      
//       // Calcola la forza
//       for(int ord = 0; ord < PTORD-1; ord++){
// 	for(int k = 0; k < 2; k++){
// 	  for(int j = 0; j < 2; j++){
	    
// 	    for(int alpha = 0; alpha < dim; alpha++){  
// 	      W.ptU[ord+1].whr[j+2*k] += ~(XiTmp.psi[alpha].whr[j])*
// 		PsiTmp.psi[alpha].ptCV[ord].whr[k];
// 	    } // alpha
	    
// 	  } // colore, j
// 	} // colore, k
//       } // ord
      
//       // Step evoluzione
//       W = W*dag(Umu.W[site_c].U[mu]);
//       W.reH();
//       //Umu.W[site_c].U[mu]= exp(act_pars.tau_f*W)*Umu.W[site_c].U[mu];
//       Wfld.W[site_c].U[mu] += act_pars.tau_f*W;

//     } // mu
 
//   } // loop sui siti

// }


  
void zero_modes_subtraction(ptBoson_fld& Umu){
  // Normalizzo il modo nullo
  for (int i1=0; i1 < dim; i1++){
    MM[i1] *= act_pars.rVol;
  }
  
  // Campo nell'algebra e sottrazione dei momenti nulli
  // nota: non mi interessa l'ordine con cui lo faccio

  #ifndef __PARALLEL_OMP__
  for(int i = 0; i < act_pars.iVol; i++){
    for(int mu = 0; mu < dim; mu++){
      W = log(Umu.W[i].U[mu]);
      W -= MM[mu];
      W.reH();
      Umu.W[i].U[mu]  = exp(W);
    }
  }
  #else

  int curr;
#pragma omp parallel private(curr, tid) num_threads(NTHR) 
  {
    tid = omp_get_thread_num();
    for(int y0 = thr_pars[tid].xi[0]; y0 < thr_pars[tid].xf[0]; y0++){
      for(int y1 = thr_pars[tid].xi[1]; y1 < thr_pars[tid].xf[1]; y1++){
	for(int y2 = thr_pars[tid].xi[2]; y2 < thr_pars[tid].xf[2]; y2++){
	  for(int y3 = thr_pars[tid].xi[3]; y3 < thr_pars[tid].xf[3]; y3++){

	    curr = y3 + act_pars.sz[3]*(y2 + act_pars.sz[2]*(y1 + act_pars.sz[1]*y0) );

	    for(int mu = 0; mu < dim; mu++){
	      Ww[tid] = log(Umu.W[curr].U[mu]);
	      Ww[tid] -= MM[mu];
	      Ww[tid].reH();
	      Umu.W[curr].U[mu]  = exp(Ww[tid]);
	    } // mu

	  } // y3
	} // y2
      } // y1
    } // y0
    
  } // parallel

#endif

}



void stochastic_gauge_fixing(ptBoson_fld& Umu){

#ifndef __PARALLEL_OMP__

  for(curr = 0; curr < act_pars.iVol; curr++){

    // preparo il sito corrente
    site_c = Umu.Z->get(curr);
    //    site_c = Umu.Z->L[curr][dim];
    W.zero();
    
    // calcolo il DmuUmu e lo metto nell'algebra
    for(int mu = 0; mu < dim; mu++){
      W += Umu.W[site_c].U[mu] - Umu.W[Umu.Z->L[curr][mu]].U[mu];
    }
    W -= dag(W);
    W.Trless();
    
    // Preparo le trasformazioni di gauge
    W1 = exp( act_pars.alpha*W);
    W2 = exp(-act_pars.alpha*W);
    
    // Moltiplico sul link corrente a sx e sul precedente a dx
    for(int mu = 0; mu < dim; mu++){
      Umu.W[site_c].U[mu] = W1*Umu.W[site_c].U[mu];
      Umu.W[Umu.Z->L[curr][mu]].U[mu] = Umu.W[Umu.Z->L[curr][mu]].U[mu]*W2;
    }
    
  } // loop sui siti

#else

#pragma omp parallel private(curr,tid,site_c)  num_threads(NTHR) 
  {
    tid = omp_get_thread_num();
    
    for(int y0 = thr_pars[tid].xi[0]; y0 < thr_pars[tid].xf[0]; y0++){
      for(int y1 = thr_pars[tid].xi[1]; y1 < thr_pars[tid].xf[1]; y1++){
	for(int y2 = thr_pars[tid].xi[2]; y2 < thr_pars[tid].xf[2]; y2++){
	  for(int y3 = thr_pars[tid].xi[3]; y3 < thr_pars[tid].xf[3]; y3++){

	    curr = y3 + act_pars.sz[3]*(y2 + act_pars.sz[2]*(y1 + act_pars.sz[1]*y0) );
	    
	    site_c = Umu.Z->L[curr][dim];
	    
	    // calcolo il DmuUmu e lo metto nell'algebra
	    Ww[tid].zero();
	    for(int mu = 0; mu < dim; mu++){
	      Ww[tid] += Umu.W[site_c].U[mu] - Umu.W[Umu.Z->L[curr][mu]].U[mu];
	    }

	    // Ww[tid] = Umu.W[site_c].U[0] - Umu.W[Umu.Z->L[curr][0]].U[0];
	    // for(int mu = 1; mu < dim; mu++){
	    //   Ww[tid] += Umu.W[site_c].U[mu] - Umu.W[Umu.Z->L[curr][mu]].U[mu];
	    // }

	    Ww[tid] -= dag(Ww[tid]);
	    Ww[tid].Trless();
	    
	    // Preparo le trasformazioni di gauge
	    Ww1[tid] = exp( act_pars.alpha*Ww[tid]);
	    Ww2[tid] = exp(-act_pars.alpha*Ww[tid]);
	    
	    // Moltiplico sul link corrente a sx e sul precedente a dx
	    for(int mu = 0; mu < dim; mu++){
	      Umu.W[site_c].U[mu] = Ww1[tid]*Umu.W[site_c].U[mu];
	      Umu.W[Umu.Z->L[curr][mu]].U[mu] = Umu.W[Umu.Z->L[curr][mu]].U[mu]*Ww2[tid];
	    }

#if ntz > 1
#pragma omp barrier
#endif
	  } // end z
	  
#if nty > 1
#pragma omp barrier
#endif
	} // end y
	
#if ntx > 1
#pragma omp barrier
#endif
      } // end x
      
#if ntt > 1
#pragma omp barrier
#endif
    }  // end t


#pragma omp barrier
  } // end parallel
#endif

}




// Beat step di evoluzione quenched
void NsptEvolve(ptBoson_fld& Umu){

#ifdef __TIMING__
  Time.reset();
#endif

  for( int t2 = 0; t2 < nspt_pars.Beat; t2++ ){
    
    // Inizializza a 0 le viariabili di appoggio, le placchette
    // e la norma
    for (int i1=0; i1 < dim; i1++){
      MM[i1].zero();
    }
    for (int i1=0; i1 <= PTORD; i1++){
      w1[i1] = 0.0;
      norm[i1]  = 0.0;
      for (int mu=0; mu < dim; mu++){
	trU[mu][i1]   = 0.0;
      }
#if GAUGE_ACTION != WIL
      w2[i1] = 0.0;
#endif
    }
  
    // Evoluzione della parte gluonica
#ifdef __TIMING__
    Time.tic_g();
#endif

    gauge_wilson(Umu);

#ifdef __TIMING__
    Time.toc_g();
#endif

//     // Sottrazione dei momenti nulli
#ifdef __TIMING__
    Time.tic_zm();
#endif

      zero_modes_subtraction(Umu);

#ifdef __TIMING__
    Time.toc_zm();
#endif

//     // Gauge fixing stocastico
#ifdef __TIMING__
    Time.tic_gf();
#endif

     stochastic_gauge_fixing(Umu);


#ifdef __TIMING__
    Time.toc_gf();
    Time.reduce();
#endif

  }

#ifdef __TIMING__
  Time.out();
#endif

}





// // Beat step di evoluzione unquenched
// void NsptEvolve(ptBoson_fld& Umu, ptSpinColor_fld& Pmu, SpinColor_fld& Xi, ptBoson_fld& Wfld){


//   for( int t2 = 0; t2 < nspt_pars.Beat; t2++ ){
    
//     // Inizializza a 0 le viariabili di appoggio, le placchette
//     // e la norma
//     for (int i1=0; i1 < dim; i1++){
//       MM[i1].zero();
//     }
//     for (int i1=0; i1 <= PTORD; i1++){
//       w1[i1] = 0.0;
//       norm[i1]  = 0.0;
// #if GAUGE_ACTION != WIL
//       w2[i1] = 0.0;
// #endif
//     }
  
//     // Forze della parte gluonica
//     gauge_wilson(Umu, Wfld);

//     // Forze della parte fermionica
//     fermion_wilson(Umu, Pmu, Xi, Wfld);
    
//     // Evoluzione
//     // Evoluzione
//     for(curr = 0; curr < act_pars.iVol; curr++){
//       site_c = Umu.Z->L[curr][dim];
//       for(int mu = 0; mu < dim; mu++){
// 	Umu.W[site_c].U[mu] = exp(Wfld.W[site_c].U[mu])*Umu.W[site_c].U[mu];
// 	// Contributo del momento nullo
// 	MM[mu] += log(Umu.W[site_c].U[mu]);		
//       }
//     }

//     // Sottrazione dei momenti nulli
//     zero_modes_subtraction(Umu);
    
//     // Gauge fixing stocastico
//     stochastic_gauge_fixing(Umu);

//   }
// }


// Misura a configurazione congelata la placchetta
void plaquette_measure(ptBoson_fld& Umu, act_params_t& act_pars){

#ifndef __PARALLEL_OMP__

  for(int i1 = 0; i1 < PTORD+1; i1++){
    w[i1] = 0.0;
  }					
  
  Cplx* app = new Cplx[PTORD+1];
  
  for(int i = 0; i < act_pars.iVol;i++){
    W1.zero();				
    link_c = Umu.Z->L[i][dim];
    for(int mu = 0; mu < dim; mu++){
      W.zero();				
      for(int nu = 0; nu < dim; nu++){
	if(nu != mu ){
	  W += Umu.staple(i, mu, nu);
	}
      }
      W1 = Umu.W[link_c].U[mu]*W;
      W1.Tr(app);
      for(int i1 = 0; i1 <= PTORD; i1++){
	w[i1] += app[i1];
      }
    }	
  }
     
  for(int i1 = 0; i1 <= PTORD; i1++){
    w[i1] /= (act_pars.iVol*48);		
  }					

  W.zero();
  W1.zero();
  delete [] app;

#else

  for(int nt = 0; nt < NTHR; nt++){
    for(int i1 = 0; i1 < PTORD+1; i1++){
      ww[nt][i1] = 0.0;
    }				
  }
  //  Cplx* app = new Cplx[PTORD+1];
  
#pragma omp parallel private(tid,link_c) num_threads(NTHR)
  {
    tid = omp_get_thread_num();

    int chunk = act_pars.iVol/NTHR;
    for(int site_x = chunk*tid; site_x < chunk*(tid+1);site_x++){
      Ww1[tid].zero();
      //      link_c = Umu.Z->L[site_x][dim];
      for(int mu = 0; mu < dim; mu++){
 	Ww[tid].zero();
 	for(int nu = 0; nu < dim; nu++){
 	  if(nu != mu ){
 	    Ww[tid] += Umu.staple(site_x, mu, nu);
 	  }
 	}
	// 	Ww1[tid] = Umu.W[link_c].U[mu]*Ww[tid];
 	Ww1[tid] = Umu.W[site_x].U[mu]*Ww[tid];
 	Ww1[tid].Tr(ww1[tid]);
	
 	for(int i1 = 0; i1 <= PTORD; i1++){
 	  ww[tid][i1] += ww1[tid][i1];
 	}
	
      }// altro link, nello stesso sito	
    } // fine siti
    
    Ww[tid].zero();
    Ww1[tid].zero();
  } // end parallel
#pragma omp barrier

  for(int i1 = 0; i1 <= PTORD; i1++){
    w[i1] = 0;
    for(int nt = 0; nt < NTHR; nt++){
      //      std::cout << ww[nt][i1].re << "\t";
      w[i1] += ww[nt][i1];
    }
    //    std::cout << std::endl;
    w[i1] /= (act_pars.iVol*72);		
  }					
  std::cout << std::endl;
 #endif

}


// Confronta 2 placchette. Se sono uguali ritorna 0, altrimenti
// la differenza tra le due
int plaquette_check(Cplx* w1, Cplx* w2){
  return ( memcmp(w1, w2, (PTORD+1)*sizeof(Cplx)) ); 
//   for(int i1 = 0; i1 <= PTORD;i1++){
//     if( (w1[i1]-w2[i1]).mod() > 1e-7 ) return 1;
//   }
//   return 0;
}


SU2 SU2rand(MyRand &Rand){
  SU2 P;
  double g;

  g = Rand.generate_gauss();
  P.whr[1].im = g;
  P.whr[2].im = g;

  g = Rand.generate_gauss();
  P.whr[1].re =  g;
  P.whr[2].re = -g;

  g = Rand.generate_gauss();
  P.whr[0].im =  g;
  P.whr[3].im = -g;

  return P;
}







#ifdef __WILLOOP_2x2__


void WL2x2(ptBoson_fld& Umu) {
  
#ifndef __PARALLEL_OMP__
  ptSU2 tmp, tmp1;
#else
  ptSU2 tmp[NTHR], tmp1[NTHR];
#endif
  
  int a,b,d,e;
  

#ifndef __PARALLEL_OMP__
  for(int site_x = 0; site_x < act_pars.iVol; site_x++) {

    for(int mu = 3; mu > 0; --mu){
      a = Umu.Z->L[site_x][5+mu];
      b = Umu.Z->L[a][5+mu];
      tmp = Umu.W[site_x].U[mu] * Umu.W[a].U[mu];

      for( int nu = 0; nu < mu; ++nu){
	d = Umu.Z->L[site_x][5+nu];
	e = Umu.Z->L[d][5+nu];
	tmp1 += tmp * ( Umu.W[b].U[nu] * Umu.W[ Umu.Z->L[b][5+nu] ].U[nu] * 
			dag( Umu.W[site_x].U[nu] * Umu.W[d].U[nu] * Umu.W[e].U[mu] * 
			     Umu.W[ Umu.Z->L[e][5+mu] ].U[mu] ));

      } // nu
    } // mu
  } // site

  for( int i1 = 1; i1 < PTORD; i1+=2){
    loop2x2 << (tmp1.ptU[i1].Tr()).re*act_pars.rVol*D12 << "\t";
  }
  loop2x2 << std::endl;

#else
#pragma omp parallel private(tid,a,b,d,e)  num_threads(NTHR) 
  {
    tid = omp_get_thread_num();
    int chunk = act_pars.iVol/NTHR;
    for(int site_x = chunk*tid; site_x < chunk*(tid+1);site_x++){
      
      for(int mu = 3; mu > 0; --mu){
	a = Umu.Z->L[site_x][5+mu];
	b = Umu.Z->L[a][5+mu];
	tmp[tid] = Umu.W[site_x].U[mu] * Umu.W[a].U[mu];
	
	for( int nu = 0; nu < mu; ++nu){
	  d = Umu.Z->L[site_x][5+nu];
	  e = Umu.Z->L[d][5+nu];
	  tmp1[tid] += tmp[tid] * ( Umu.W[b].U[nu] * Umu.W[ Umu.Z->L[b][5+nu] ].U[nu] * 
				    dag( Umu.W[site_x].U[nu] * Umu.W[d].U[nu] * Umu.W[e].U[mu] * 
					 Umu.W[ Umu.Z->L[e][5+mu] ].U[mu] ));
	} // nu
      } // mu
      
    } // end site_x
    
  } // end parallel
  

  // Reduce
  for(int tid = 1; tid < NTHR; tid++){
    tmp1[0] += tmp1[tid];
  }

  // std::cout << std::endl;
  for( int i1 = 1; i1 <= PTORD; i1+=2){
    // std::cout << (tmp1[0].ptU[i1].Tr()).re*act_pars.rVol*D18 << "\t";
    loop2x2 << (tmp1[0].ptU[i1].Tr()).re*act_pars.rVol*D12 << "\t";
  }

  // std::cout << std::endl;
  loop2x2 << std::endl;

#endif

}

#endif
