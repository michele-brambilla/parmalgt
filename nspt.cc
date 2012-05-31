#include "nspt.h"
#include "GammaStuff.h"
#include <iomanip>

ptSU3 W, W1, W2, Wsym;
ptSU3 *MM = new ptSU3[dim];

int link_c, site_c, site_up, xx[4], curr;
const int offset = dim+1;

extern Cplx* w;  // usero' come variabile d'appoggio
extern Cplx* w1; // usero' per placchetta 1x1
extern Cplx* w2; // usero' per placchetta 2x1
extern double* norm;
extern double* norm1;
extern Cplx** trU;
extern int ratio; 
extern ptSU3* U;
extern nspt_params_t nspt_pars;
extern act_params_t  act_pars;
extern thread_params_t thr_pars[NTHR];

extern std::ofstream plaqfile, normfile, logfile, trufile;

#ifdef __WILLOOP_2x2__
extern std::ofstream loop2x2;
#endif


extern MyRand Rand[NTHR];

#ifdef __PARALLEL_OMP__
ptSU3 Ww[NTHR], Ww1[NTHR], Ww2[NTHR], Wwsym[NTHR];
ptSU3 MMm[NTHR][dim];
SU3 P[NTHR];
Cplx ww[NTHR][ORD+1], ww1[NTHR][ORD+1], ww2[NTHR][ORD+1];
double normm[NTHR][ORD+1];
Cplx ttrU[NTHR][dim][ORD];
int tid;

#endif

// for (single step of) Fourier Accelerated gauge fixing
int coord[dim];
ptSU3_fld* Wgauge;
fftw_plan *planFA;

// for momentum space analysis of the norm
fftw_plan *planGluon;

#ifdef __K_MOM_ANALYSIS__
double *knorm[ORD];
extern FILE *FPkn;
#endif

#ifdef __TIMING__
extern PRlgtTime Time;
#endif

void gauge_wilson(ptGluon_fld& Umu){

#ifndef __PARALLEL_OMP__

  for(curr = 0; curr < act_pars.iVol; curr++){

#ifdef SFBC
    int t = curr / (act_pars.sz[1]*act_pars.sz[2]*act_pars.sz[3]); 
#endif
    
    for(int mu = 0; mu < dim; mu++){

      // DH, Feb. 6, 2012
      // FIXME: THIS IS JUST A ROUGH HACK!
      // 
      // The boundary conditions in the LNWW paper (hep-lat/9207009)
      // are:
      //       U_k(x) = V_k(\vec x),    @ x_0 = 0
      //       U_k(x) = V'_k('vec x),   @ x_0 = T
      //
      // This means, that we can use the following hack to get us some
      // nice SF boundary conditions. We simply assume that
      // act_pars.sz[0] == (T+1) instead of T. Thus, we may
      // 1) We initialize U to V(t), getting the above automatically
      //    right
      // 2) in the update for U_mu(x, t):
      //        if we have t == T or (t == 0 and mu != 0)  
      //          ==> DO NOTHING!
#ifdef SFBC
      if ( (!t && mu) || (t == act_pars.sz[0] - 1) ) continue;
#endif

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
      for (int i1=1; i1 <= ORD; i1++){		
	w1[i1] += w[i1];
	norm[i1-1] += ((U[link_c][i1-1]*dag(U[link_c][i1-1])).Tr()).re;
	trU[mu][i1-1] += U[link_c][i1-1].Tr();
      }

      
      // Se ho azioni improved faccio lo stesso sulle
      // parti di improvement
#if GAUGE_ACTION != WIL
      W2 = U[link_c]*W2;			
      W2.Tr(w);
      for (int i1=0; i1 <= ORD; i1++){		
	w2[i1] += w[i1];
      }
      
      W1 = act_pars.c0*W1 + act_pars.c1*W2;
#endif

      // DH,  Feb 6, 2012
      // The original code looked like this

      // Costringo nell'algebra
      //W1 = W1.reH();					
      //W1 *= act_pars.tau_g;
      
      // Aggiungo fluttuazione al primo ordine
      //W1[0] -= act_pars.stau*SU3rand(Rand[0]);
      //U[link_c] = exp(W1)*U[link_c];	
      
      // Contributo del momento nullo
      //MM[mu] += log(U[link_c]);		
      
      // Now we have this

      // DISCLAIMER
      //
      // AS OF TODAY, I DID NOT TEST THE NON-PARALLEL VERSION OF THE CODE
      // USE AT YOUR OWN RISK (AND PLEASE REMOVE THIS NOTE IF YOU TRY AND 
      // SEE THAT IT WORKS, THANKS!).

      ptsu3 tmp = W1.reH() * act_pars.tau_g; // take W1 to the algebra
      tmp[0]  -= act_pars.stau*SU3rand(Rand[0]); // add noize
      U[link_c] = exp<bgf::AbelianBgf>(W1)*U[link_c]; // back to the group
      MM[mu] += get_q(U[link_c]); //
      
    } //end mu
    
  } // fine siti


#else
  
  for(tid = 0;tid < NTHR; tid++) {
    for(int ord = 0; ord < ORD; ord++){
      ww[tid][ord]  = 0;
      ww1[tid][ord] = 0;
      ww2[tid][ord] = 0;
      normm[tid][ord] = 0;
      for( int mu = 0; mu < dim; mu++)
	{
	  ttrU[tid][mu][ord] = Cplx(0,0);
	}
    }
    ww[tid][ORD]  = 0;
    ww1[tid][ORD] = 0;
    ww2[tid][ORD] = 0;
    
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
              // DH, Feb. 6, 2012, c.f. comment above

              // SCHRÖDINGER FUNCTIONAL BOUNDARY CONDITIONS
              // PROCEED WITH CARE!!!

              // READ THE COMMENT IN THE SINGLE THREADED VERSION ABOVE!
#ifdef SFBC
              // y0== t!
              if ( (!y0 && mu) || (y0 == act_pars.sz[0] - 1) ) continue;
              // no seriously, read the comment!!
#endif
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
	      for (int i1=1; i1 <= ORD; i1++){
		ww1[tid][i1] += ww[tid][i1];

		normm[tid][i1-1] += ((U[link_c][i1-1]*dag(U[link_c][i1-1])).Tr()).re;

		ttrU[tid][mu][i1-1] += U[link_c][i1-1].Tr();
	      }

	      // Se ho azioni improved faccio lo stesso sulle
	      // parti di improvement
#if GAUGE_ACTION != WIL
	      Ww2[tid] = U[link_c]*Ww2[tid];			
	      Ww2[tid].Tr(ww[tid]);
	      for (int i1=0; i1 <= ORD; i1++){		
		ww2[tid][i1] += ww[tid][i1];
	      }
      
	      Ww1[tid] = Ww1[tid]*act_pars.c0 + Ww2[tid] * act_pars.c1;
#endif
              
              // DH Feb. 6, 2012
              ptsu3 tmp  = Ww1[tid].reH(); // take to the algebra
              tmp *= act_pars.tau_g; // multipy by tau_g
              tmp[0] -= act_pars.stau*SU3rand(Rand[tid]); // add noise
              U[link_c] = exp<bgf::AbelianBgf, ORD>(tmp)*U[link_c]; // back to SU3
              MMm[tid][mu] += get_q(U[link_c]); // zero momentum contribution
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
  
    for(int ord = 0; ord < ORD; ord++){
      w1[ord] += ww1[tid][ord];
      w2[ord] += ww2[tid][ord];
      norm[ord] += normm[tid][ord];
      for(int mu = 0; mu < dim; mu++){
	trU[mu][ord] += ttrU[tid][mu][ord];
      }
    }
    w1[ORD] += ww1[tid][ORD];
    w2[ORD] += ww2[tid][ORD];

    for(int mu = 0; mu < dim; mu++){
      MM[mu] += MMm[tid][mu];
    }
  }

#endif

  
  for (int i1=0; i1 < ORD; i1++){		
#if GAUGE_ACTION == WIL
    plaqfile << w1[i1].re/(double)(72*act_pars.iVol) << "\t" 
	     << w1[i1].im/(double)(72*act_pars.iVol) << std::endl;
#else
    plaqfile << w1[i1].re/(double)(72*act_pars.iVol) << "\t" 
	     << w1[i1].im/(double)(72*act_pars.iVol) << "\t"
	     << w2[i1].re/(double)(72*act_pars.iVol) << "\t" 
	     << w2[i1].im/(double)(72*act_pars.iVol) << std::endl;
#endif
    normfile << norm[i1]/(double)(3*dim*act_pars.iVol) << std::endl;
    for( int mu = 0; mu < dim; mu++)
      {
	trufile  << trU[mu][i1].re /(double)(3*act_pars.iVol) << "\t" 
		 << trU[mu][i1].im /(double)(3*act_pars.iVol) << "\t";
      }
    trufile  << std::endl;
  }

#if GAUGE_ACTION == WIL
  plaqfile << w1[ORD].re/(double)(72*act_pars.iVol) << "\t" 
	   << w1[ORD].im/(double)(72*act_pars.iVol) << std::endl;
#else
  plaqfile << w1[ORD].re/(double)(72*act_pars.iVol) << "\t" 
	   << w1[ORD].im/(double)(72*act_pars.iVol) << "\t"
	   << w2[ORD].re/(double)(72*act_pars.iVol) << "\t" 
	   << w2[ORD].im/(double)(72*act_pars.iVol) << std::endl;
#endif

} // fine evoluzione gluoni





SpinColor ZeroMom;
SpinColor XiTmp;
ptSpinColor PsiTmp;
#ifdef __PARALLEL_OMP__
SpinColor ZM[NTHR];
SpinColor XT[NTHR];
ptSpinColor PT[NTHR];
#endif

void fermion_wilson(ptGluon_fld& Umu, ptSpinColor_fld& Pmu, SpinColor_fld& Xi){

  memset((void*)&ZeroMom, 0, sizeof(SpinColor));

#ifndef __PARALLEL_OMP__
  // genero campo gaussiano

  Xi.gauss(Rand[0]);

  // Sottraggo il momento nullo
  for( int iv = 0; iv < act_pars.iVol; iv++){
    ZeroMom += Xi[iv];
  }
  ZeroMom *= act_pars.rVol;

  for( int iv = 0; iv < act_pars.iVol; iv++){
    Xi[iv] -= ZeroMom;
    Pmu.scfld->psi[iv] = Xi[iv];
  }

#else

  // Sottraggo il momento nullo
  memset((void*)ZM, 0, NTHR*sizeof(SpinColor));
   
#pragma omp parallel private(ZeroMom) num_threads(NTHR)
  { 
    int tid = omp_get_thread_num();
    Xi.gauss(Rand,tid);
#pragma omp for
    for(int iv = 0; iv < act_pars.iVol; iv++){
      ZM[tid] += Xi[iv];
    }
  }
  
  ZeroMom = ZM[0];
  for(int nt = 1; nt < NTHR; nt++){
    ZeroMom += ZM[nt];
  }
  ZeroMom *= act_pars.rVol;
  
#pragma omp parallel for num_threads(NTHR)
  for( int iv = 0; iv < act_pars.iVol; iv++){
    Xi[iv] -= ZeroMom;
    Pmu.scfld->psi[iv] = Xi[iv];
  } // pragma omp parallel for -> now single thread
  

#endif

  memset((void*)Pmu.psi, 0, act_pars.iVol*sizeof(ptSpinColor));
  Pmu.fillPT(ORD-2,Umu);

#ifndef __PARALLEL_OMP__
  // cicla su tutti i siti
  for(curr = 0; curr < act_pars.iVol; curr++){
    site_c = Umu.Z->L[curr][dim];

    for( int mu = 0; mu < dim; mu++ ){
      
      // variabile d'appoggio
      W.zero();
      XiTmp   = Xi[Umu.Z->L[curr][offset+mu]];
      // PsiTmp  = Pmu[site_c];
      // PsiTmp += PsiTmp.gmleft(mu);

      for(int nu = 0; nu < dim; nu++)
	{
	  PsiTmp[nu] = Pmu[site_c][nu] + gmuval[mu][nu]*Pmu[site_c][gmuind[mu][nu]];
	}

      // Calcola la forza
      for(int ord = 0; ord < ORD-1; ord++){
	for(int k = 0; k < 3; k++){
	  for(int j = 0; j < 3; j++){
	    
	    for(int alpha = 0; alpha < dim; alpha++){  
	      W[ord+1].whr[j+3*k] += ~(XiTmp[alpha].whr[j])*
		PsiTmp[alpha][ord].whr[k];
	    } // alpha
	    
	  } // colore, j
	} // colore, k
      } // ord
      
      // Step evoluzione
      W = W*dag(Umu.W[site_c][mu]);
      Umu.W[site_c][mu]= exp<bgf::AbelianBgf,ORD>(act_pars.tau_f*W.reH())*Umu.W[site_c][mu];
    } // mu

  } // loop sui siti
#else
#pragma omp parallel private(tid,site_c,curr) num_threads(NTHR)
  {
    tid = omp_get_thread_num();

    for(int y0 = thr_pars[tid].xi[0]; y0 < thr_pars[tid].xf[0]; y0++){
      for(int y1 = thr_pars[tid].xi[1]; y1 < thr_pars[tid].xf[1]; y1++){
	for(int y2 = thr_pars[tid].xi[2]; y2 < thr_pars[tid].xf[2]; y2++){
	  for(int y3 = thr_pars[tid].xi[3]; y3 < thr_pars[tid].xf[3]; y3++){

	    curr = y3 + act_pars.sz[3]*(y2 + act_pars.sz[2]*(y1 + act_pars.sz[1]*y0) );
	    
	    site_c = Umu.Z->L[curr][dim];
	    
	    for( int mu = 0; mu < dim; mu++ ){
	      
	      // variabile d'appoggio
	      Ww[tid].zero();
	      XT[tid]   = Xi[Umu.Z->L[curr][offset+mu]];
	      // PT[tid]  = Pmu[site_c];
	      // PT[tid] += PT[tid].gmleft(mu);

	      for(int nu = 0; nu < dim; nu++)
		{
		  PT[tid][nu] = Pmu[site_c][nu] + 
                    Pmu[site_c][gmuind[mu][nu]] * gmuval[mu][nu];
		}
	      
	      // Calcola la forza
	      for(int ord = 0; ord < ORD-1; ord++){
		for(int k = 0; k < 3; k++){
		  for(int j = 0; j < 3; j++){
		    
		    for(int alpha = 0; alpha < dim; alpha++){  
		      Ww[tid][ord+1].whr[j+3*k] += ~(XT[tid][alpha].whr[j])*
			PT[tid][alpha][ord].whr[k];
		    } // alpha
		    
		  } // colore, j
		} // colore, k
	      } // ord
	      
	      // Step evoluzione
	      Ww[tid] = Ww[tid]*dag(Umu.W[site_c][mu]);
	      Umu.W[site_c][mu]= exp<bgf::AbelianBgf,ORD>(Ww[tid].reH() * act_pars.tau_f)
                *Umu.W[site_c][mu];
	    } // mu
	    
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
  } // parallel
#endif

}


  
void zero_modes_subtraction(ptGluon_fld& Umu){
  // Normalizzo il modo nullo
  for (int i1=0; i1 < dim; i1++){
    MM[i1] *= act_pars.rVol;
  }
  
  
  // Campo nell'algebra e sottrazione dei momenti nulli
  // nota: non mi interessa l'ordine con cui lo faccio

  #ifndef __PARALLEL_OMP__
  for(int i = 0; i < act_pars.iVol; i++){
    for(int mu = 0; mu < dim; mu++){
      W.ptU() = log(Umu.W[i][mu]);
      W -= get_q(MM[mu]);
      W.ptU() = W.reH();
      W.bgf().reH();
      Umu.W[i][mu]  = exp<bgf::AbelianBgf,ORD>(W.ptU());
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


#ifdef SFBC // EXPERIMENTAL SF BOUNDARY STUFF
            // if we proceed this way, we omit the gauge fixing for q_0(x)_{x_0  = 0}
            // alltogether, c.f. the note below.
            // However, this seems to work
            if ( !y0 || (y0 == act_pars.sz[0] - 1) ) continue;
#endif

	    for(int mu = 0; mu < dim; mu++){
              // DH Feb. 1, 2012
              // TODO: This is more a dirty hack than anything else.
              // The problem is that the old log function silently
              // assumed that we have U = 1 + ..., which is not true
              // anymore. Thus, we have an extra constructor call
              // here... This should be fixed!
              //Umu.W[curr][mu].bgf() = bgf::unit();
	      //Ww[tid] = log(Umu.W[curr][mu]);
	      //Ww[tid] -= MM[mu];
	      //Ww[tid].reH();
	      //Umu.W[curr][mu]  = exp(Ww[tid]);


              //Ww[tid].ptU() = get_q(Umu.W[curr][mu]);
              //Ww[tid] = -get_q(MM[mu]);
              
              //Ww[tid].bgf().reH();
              //Umu.W[curr][mu].ptU() = exp<bgf::AbelianBgf, ORD>(Ww[tid].ptU()).ptU();
              //for (int i = 0; i < ORD; ++i)
              //  Umu.W[curr][mu][i] = Umu.W[curr][mu].bgf().ApplyFromRight(Umu.W[curr][mu][i]);
              //Umu.W[curr][mu] *= exp<bgf::AbelianBgf,
              //ORD>(MM[mu].ptU()*-1);
              Umu.W[curr][mu] = exp<bgf::AbelianBgf, ORD>(-1*MM[mu].reH())*Umu.W[curr][mu];
	    } // mu

	  } // y3
	} // y2
      } // y1
    } // y0
    
  } // parallel

#endif

}








void stochastic_gauge_fixing(ptGluon_fld& Umu){

#ifndef __PARALLEL_OMP__
  for(curr = 0; curr < act_pars.iVol; curr++){

    // preparo il sito corrente
    site_c = Umu.Z->get(curr);
    //    site_c = Umu.Z->L[curr][dim];
    W.zero();
    
    // calcolo il DmuUmu e lo metto nell'algebra
    for(int mu = 0; mu < dim; mu++){
      W += Umu.W[site_c][mu] - Umu.W[Umu.Z->L[curr][mu]][mu];
    }

    // DH Feb. 6 2012
    
    // Apart from the original formula (in //// comments), whe had the
    // instructions W -= dag(w) and W.Trless(). This roughly amounts
    // to reH() and I used this below!
    //
    // DISCLAIMER: I HAVE NOT TRIED THE NON-PARALLEL VERSION
    //             CONSIDER IT TO BE BROKEN UNLESS YOU CHECK IT!!!

    //W -= dag(W);
    
    //// Forma alternativa
    //// for(int mu = 0; mu < dim; mu++){
    ////   W += log(Umu.W[site_c][mu])- log(Umu.W[Umu.Z->L[curr][mu]][mu]);
    //// }

    //W.Trless();
    
    // Preparo le trasformazioni di gauge

    // The new exp takes arguments in the algebra, thus we use reH().

    W1 = exp<bgf::AbelianBgf, ORD>( act_pars.alpha*W.reH());
    W2 = exp<bgf::AbelianBgf, ORD>(-act_pars.alpha*W.reH());


    // Moltiplico sul link corrente a sx e sul precedente a dx
    for(int mu = 0; mu < dim; mu++){
      Umu.W[site_c][mu] = W1*Umu.W[site_c][mu];
      Umu.W[Umu.Z->L[curr][mu]][mu] = Umu.W[Umu.Z->L[curr][mu]][mu]*W2;
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
	    // DH, Feb. 7 2012
#ifdef SFBC // EXPERIMENTAL SF BOUNDARY STUFF
            // if we proceed this way, we omit the gauge fixing for q_0(x)_{x_0  = 0}
            // alltogether, c.f. the note below.
            // However, this seems to work
            if ( !y0 || (y0 == act_pars.sz[0] - 1) ) continue;
#endif
    
	    // calcolo il DmuUmu e lo metto nell'algebra
	    Ww[tid].zero();
            
            /*/
            for(int mu = 0; mu < 4; mu++)
              Ww[tid] += get_q(Umu.W[site_c][mu]) - get_q(Umu.W[Umu.Z->L[curr][mu]][mu]);
            /*/
            for(int mu = 0; mu < 4; mu++){
              // Calculate
              // U_mu(x) * V^\dagger_mu(x - mu) 
              // * U_mu ^\dagger (x - mu) * V_mu(x - mu)
              ptSU3 Udag = dag(Umu.W[Umu.Z->L[curr][mu]][mu]);
              ptSU3 &U = Umu.W[site_c][mu];
              bgf::AbelianBgf Vdag = U.bgf().dag(), V = Udag.bgf().dag();
              Ww[tid] += U*Vdag*Udag*V;
              //Ww[tid].bgf() = bgf::AbelianBgf();
              //ptSU3 one = Ww[tid]*inverse(Ww[tid]);
              //
              //std::cout << "xxxxxxxxxxxxxxxxxxx" << one.bgf() <<
              //std::endl;
              //for (int s = 0; s < ORD; ++s)
              //  std::cout << one[s].Norm() << std::endl;
              //
              //for (int nu = 0; nu < 4; ++nu){
              /*/

                std::vector<std::vector<Cplx> > plaq_tmp(ORD, std::vector<Cplx>(16));
                for (int mu2 = 0; mu2 < 4; ++mu2)
                  for (int nu2 = 0; nu2 < 4; ++nu2)
                    for (int i = 0; i < ORD; ++i)
                      plaq_tmp[i][mu2*4 + nu2] =
	    (Umu.W[site_c][mu2]*Umu.staple(curr, mu2, nu2))[i].Tr();
            /*/
              // Umu.W[site_c][nu] = Ww[tid]*Umu.W[site_c][nu];
              // Umu.W[Umu.Z->L[curr][nu]][nu] =
              //   Umu.W[Umu.Z->L[curr][nu]][nu] * dag(Ww[tid]);
                /*/
                for (int mu2 = 0; mu2 < 4; ++mu2)
                  for (int nu2 = 0; nu2 < 4; ++nu2)
                    for (int i = 0; i < ORD; ++i){
                      Cplx p_new = (Umu.W[site_c][mu2]*Umu.staple(curr, mu2, nu2))[i].Tr();
                      Cplx p_old = plaq_tmp[i][mu2*4 + nu2];
                      Cplx tmp = p_new - p_old;
                      double eps = sqrt(tmp.re*tmp.re + tmp.im*tmp.im);
                      if (eps > 1e-14)
                        std::cout << i << ": " <<  eps << std::endl
                          ;//<< "   value: " << p_old << std::endl;
                    }
                    /*/
                //}
            }
            //*/
	    //for(int mu = 0; mu < dim; mu++){
            //  Ww[tid] += get_q(Umu.W[site_c][mu]);
            //  ptSU3 Uref = Umu.W[Umu.Z->L[curr][mu]][mu];
            //  ptsu3 q = get_q(Uref);
            //  bgf::AbelianBgf V = Uref.bgf();
            //  bgf::AbelianBgf Vinv = V.inverse();
            //  for (int i = 0; i < ORD; ++i)
            //    Ww[tid][i] -= Vinv.ApplyFromLeft(V.ApplyFromRight(q[i]));
	    //}
            //for (int i = 0; i < ORD; ++i)
            //  Ww[tid][i] += Umu.W[site_c][3][i]
            //    - Umu.W[site_c][3].bgf().ApplyFromLeft
            //    ( Umu.W[Umu.Z->L[curr][3]][3].bgf().inverse().ApplyFromRight
            //      ( Umu.W[Umu.Z->L[curr][3]][3][i]));

	    // Ww[tid] = Umu.W[site_c][0] - Umu.W[Umu.Z->L[curr][0]][0];
	    // for(int mu = 1; mu < dim; mu++){
	    //   Ww[tid] += Umu.W[site_c][mu] - Umu.W[Umu.Z->L[curr][mu]][mu];
	    // }

            // DH Feb. 6, 2012
            // I replaced the instructions below with reH(), c.f. my
            // comment above (except for the disclaimer, I tried this out!)

	    //Ww[tid] -= dag(Ww[tid]);
	    //Ww[tid].Trless();
	    
	    // Preparo le trasformazioni di gauge
	    Ww1[tid] = exp<bgf::AbelianBgf, ORD>( Ww[tid].reH() * act_pars.alpha);
	    Ww2[tid] = exp<bgf::AbelianBgf, ORD>( Ww[tid].reH() * -act_pars.alpha);
            //Ww[tid] *= .25;
            //Ww[tid].bgf() = bgf::AbelianBgf();
            //std::cout << Ww[tid].bgf() - bgf::unit<bgf::AbelianBgf>() << std::endl;
            //Ww1[tid] = Ww[tid];// * alpha;
            //Ww1[tid].bgf() = bgf::AbelianBgf();
            //Ww2[tid] = inverse(Ww[tid]);// * alpha;
            // BEWARE that the above may throw an excpetion if reH
            // does not know what to do!
            /*/
            std::vector<std::vector<Cplx> > plaq_tmp(ORD, std::vector<Cplx>(16));
            for (int mu = 0; mu < 4; ++mu)
              for (int nu = 0; nu < 4; ++nu)
                for (int i = 0; i < ORD; ++i)
                  plaq_tmp[i][mu*4 + nu] = (Umu.W[site_c][mu]*Umu.staple(curr, mu, nu))[i].Tr();
            /*/
	    // Moltiplico sul link corrente a sx e sul precedente a dx
	    for(int mu = 0; mu < 4; mu++){ 
	      Umu.W[site_c][mu] = Ww1[tid]*Umu.W[site_c][mu];
	      Umu.W[Umu.Z->L[curr][mu]][mu] = Umu.W[Umu.Z->L[curr][mu]][mu]*Ww2[tid];
	    }
            /*/
            for (int mu = 0; mu < 4; ++mu)
              for (int nu = 0; nu < 4; ++nu)
                for (int i = 0; i < ORD; ++i){
                  Cplx p_new = (Umu.W[site_c][mu]*Umu.staple(curr, mu, nu))[i].Tr();
                  Cplx p_old = plaq_tmp[i][mu*4 + nu];
                  Cplx tmp = p_new - p_old;
                  double eps = sqrt(tmp.re*tmp.re + tmp.im*tmp.im);
                  if (eps > 1e-14)
                    std::cout << i << ": " <<  eps << std::endl
                      ;//<< "   value: " << p_old << std::endl;
                }
                /*/

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


#ifdef __K_MOM_ANALYSIS__

  fftw_execute(planGluon[0]);
  for(curr = 0; curr < act_pars.iVol; curr++){
    for (int i1=0; i1 < ORD; i1++){		
      knorm[i1][curr] = 0;
    }
    for(int mu = 0; mu < dim; mu++){
      for (int i1=0; i1 < ORD; i1++){		
	knorm[i1][curr] += ( Umu.W[curr][mu][i1]*
			     dag(Umu.W[curr][mu][i1]) 
			     ).Tr().re;
      }
    }
  }

  fftw_execute(planGluon[1]);
  
  for(curr = 0; curr < act_pars.iVol; curr++){
    for(int mu = 0; mu < dim; mu++){
      Umu.W[curr][mu] *= act_pars.rVol;
    }
  }

  fwrite(knorm[0],sizeof(double),act_pars.iVol*ORD ,FPkn);

#endif


}






void FAstochastic_gauge_fixing(ptGluon_fld& Umu){

#ifdef  __PARALLEL_OMP__
#pragma omp for schedule(static)
#endif
  for(int x = 0; x < act_pars.iVol; x++){
    Wgauge->W[Umu.Z->get(x)].zero();
    for(int mu = 0; mu < dim; mu++){
      Wgauge->W[Umu.Z->get(x)] += ( log(U[get(&Umu,x,mu)]) - 
      log(U[get(&Umu,x,-1,mu,mu)]) );
    }
  }
  
  fftw_execute(planFA[0]);
  
#ifdef PBC
  Wgauge->W[0] *= 0.0;
#endif  
  
#ifdef  __PARALLEL_OMP__
#pragma omp for schedule(static)
#endif
  for(int x = 1; x < act_pars.iVol; x++){
    Umu.Z->get(x,coord);
    Wgauge->W[Umu.Z->get(x)] *= -(act_pars.alpha/(Umu.Z->p2hat[0][coord[0]]+
						  Umu.Z->p2hat[1][coord[1]]+
						  Umu.Z->p2hat[2][coord[2]]+
						  Umu.Z->p2hat[3][coord[3]]));
  }
  
  
  fftw_execute(planFA[1]);
  
#ifdef  __PARALLEL_OMP__
#pragma omp for schedule(static)
#endif
  for(int x = 0; x < act_pars.iVol; x++){
    Wgauge->W[x] /= (double)act_pars.iVol;
    Wgauge->W[x] = exp<bgf::AbelianBgf, ORD>(Wgauge->W[x].reH());
  }
    
  
#ifdef  __PARALLEL_OMP__
#pragma omp for schedule(static)
#endif
  for(int x = 0; x < act_pars.iVol; x++){
    for(int mu = 0; mu < dim; mu++){
      U[get(&Umu,x,mu)] = Wgauge->get(x)*Umu.get(x,mu)*dag(Wgauge->get(x,1,mu));
    }
  }
  

#ifdef __K_MOM_ANALYSIS__

  fftw_execute(planGluon[0]);

  for(curr = 0; curr < act_pars.iVol; curr++){
    for (int i1=0; i1 < ORD; i1++){		
      knorm[i1][curr] = 0;
    }
    for(int mu = 0; mu < dim; mu++){
      for (int i1=0; i1 < ORD; i1++){		
	knorm[i1][curr] += ( Umu.W[curr][mu][i1]*
			     dag(Umu.W[curr][mu][i1]) 
			     ).Tr().re;
      }
    }
  }

  fftw_execute(planGluon[1]);
  
  for(curr = 0; curr < act_pars.iVol; curr++){
    for(int mu = 0; mu < dim; mu++){
      Umu.W[curr][mu] *= act_pars.rVol;
    }
  }

  fwrite(knorm[0],sizeof(double),act_pars.iVol*ORD ,FPkn);

#endif
  
}


//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
///
///  u0_print
///
///  Write the norm of the temporal component of U at the t = 0
///  boundary to a file. This is used to check if the naive gauge
///  fixing does something funny to this particular component.
///
///  \author Dirk Hesse <herr.dirk.hesse@gmail.com>
///  \date Wed Feb  8 15:01:36 2012

inline void u0_print(ptGluon_fld &Umu) {
  std::ofstream of("q0.nor", std::ios_base::app);
  std::vector<double> nor(ORD);
  for(int x = 0; x < act_pars.sz[0]; ++x)
    for(int y = 0; y < act_pars.sz[1]; ++y)
      for(int z = 0; z < act_pars.sz[2]; ++z){
        int curr = act_pars.sz[3]*(z + act_pars.sz[2]*(y + act_pars.sz[1]*x) );
        long link_c = get(&Umu, curr, 0);
        for (int i = 0; i < ORD; ++i)
          nor[i] += U[link_c][i].Norm();
      }
  for (int i = 0; i < ORD; ++i)
    of << nor[i]/act_pars.sz[0]/act_pars.sz[0]/act_pars.sz[0] << " ";
  of << std::endl;
  of.flush();
  of.close();
}

struct ThreadInfo {
  ThreadInfo (const thread_params_t * const tp, const int& n) :
    t (tp[n].xi[0]),
    x (tp[n].xi[1]),
    y (tp[n].xi[2]),
    z (tp[n].xi[3]),
    dt (tp[n].xf[0] - tp[n].xi[0]),
    dx (tp[n].xf[1] - tp[n].xi[1]),
    dy (tp[n].xf[2] - tp[n].xi[2]),
    dz (tp[n].xf[3] - tp[n].xi[3]) { }
  int t,x,y,z,dt,dx,dy,dz;
};

template <class C>
inline void meas_on_timeslice(const ptGluon_fld& U, const int& t, C& f){
  #pragma omp parallel num_threads(NTHR)
  {
    ThreadInfo ti(thr_pars, omp_get_thread_num());
    pt::Point n = U.mk_point(t, ti.x, ti.y, ti.z);
    for (int n1_ = 0; n1_ < ti.dx; ++n1_, n += pt::Direction::x)
      for (int n2_ = 0; n2_ < ti.dy; ++n2_, n += pt::Direction::y)
        for (int n3_ = 0; n3_ < ti.dz; ++n3_, n += pt::Direction::z)
          f(U, n);
  }
}

template <class C>
inline void apply_on_timeslice(ptGluon_fld& U, const int& t, C& f){
  #pragma omp parallel num_threads(NTHR)
  {
    ThreadInfo ti(thr_pars, omp_get_thread_num());
    pt::Point n = U.mk_point(t, ti.x, ti.y, ti.z);
    for (int n1_ = 0; n1_ < ti.dx; ++n1_, n += pt::Direction::x)
      for (int n2_ = 0; n2_ < ti.dy; ++n2_, n += pt::Direction::y)
        for (int n3_ = 0; n3_ < ti.dz; ++n3_, n += pt::Direction::z)
          f(U, n);
  }
}

// This is the plaquette at t = 0, arranged such that the derivative
// w.r.t. eta may be inserted at the very end.

struct P_lower {
  ptSU3 val;
  P_lower () : val(bgf::zero<bgf::AbelianBgf>()) { }
  void operator()(const ptGluon_fld& U, const pt::Point n){
    pt::Direction t = pt::Direction::t;
    for (pt::Direction k(1); k.good(); ++k)
      val += U(n, k).bgf() * U(n + k, t) * 
        dag( U(n +t, k) ) * dag( U(n, t) );
  }
};

// This is the plaquette at t = T, arranged such that the derivative
// w.r.t. eta may be inserted at the very end.

struct P_upper {
  ptSU3 val;
  P_upper () : val(bgf::zero<bgf::AbelianBgf>()) { }
  void operator()(const ptGluon_fld& U, const pt::Point n){
    pt::Direction t = pt::Direction::t;
    ptSU3 tmp;
    for (pt::Direction k(1); k.good(); ++k)
      val += dag( U (n + t, k).bgf() ) * dag( U(n, t) ) *
              U(n, k) * U(n + k, t);
  }
};

struct P_spatial {
  ptSU3 val;
  P_spatial () : val(bgf::zero<bgf::AbelianBgf>()) { }
  void operator()(const ptGluon_fld& U, const pt::Point n){
    for (pt::Direction k(1); k.good(); ++k)
      for (pt::Direction l(k + 1); l.good(); ++l)
        val += U(n,k) * U(n + k, l)
          * dag( U(n + l, k) ) * dag( U(n, l) );
  }
};

inline void to_bin_file(std::ofstream& of, const Cplx& c){
  of.write(reinterpret_cast<const char*>(&c.re), sizeof(double));
  of.write(reinterpret_cast<const char*>(&c.im), sizeof(double));
}

inline void write_file(const ptSU3& U, const Cplx& tree, const std::string& fname){
  std::ofstream of(fname.c_str(), std::ios_base::app |  std::ios_base::binary);
  to_bin_file(of, tree);
  for (int i = 0; i < ORD; ++i)
    to_bin_file(of, U[i].Tr());
  of.close();
}

inline void measure(const ptGluon_fld &U){
  std::vector<double> nor(ORD*2*3 + 2);
  SU3 dC; // derivative of C w.r.t. eta
  dC(0,0) = Cplx(0, -2./act_pars.sz[1]);
  dC(1,1) = Cplx(0, 1./act_pars.sz[1]);
  dC(2,2) = Cplx(0, 1./act_pars.sz[1]);
  SU3 d_dC; // derivative of C w.r.t. eta and nu
  d_dC(0,0) = Cplx(0, 0);
  d_dC(1,1) = Cplx(0, 2./act_pars.sz[1]);
  d_dC(2,2) = Cplx(0, -2./act_pars.sz[1]);

  // shorthand for T, L, V3
  int T = act_pars.sz[0] - 1;
  int L = act_pars.sz[1];
  int V3 = L*L*L;

  P_lower Pl; // Plaquettes U_{0,k} at t = 0
  P_upper Pu; // Plaquettes U_{0,k} at t = T - 1
  meas_on_timeslice(U, 0, Pl);
  meas_on_timeslice(U, T - 1, Pu);
  P_upper Pm; // Plaquettes U_{0,k} at t = 1, ..., T - 2
  for (int t = 1; t < T - 1; ++t)
    meas_on_timeslice(U, t, Pm);

  P_spatial Ps; // Spatial plaq. at t = 1, ..., T-1
  for (int t = 1; t < T; ++t)
    meas_on_timeslice(U, t, Ps);

  // Evaluate Gamma'
  ptSU3 tmp = Pl.val + Pu.val;
  for_each(tmp.begin(), tmp.end(), pta::mul(dC));
  Cplx tree = tmp.bgf().ApplyFromLeft(dC).Tr();
  write_file(tmp, tree, "Gp.bindat");

  // Evaluate v
  tmp = Pl.val + Pu.val;
  for_each(tmp.begin(), tmp.end(), pta::mul(d_dC));
  tree = tmp.bgf().ApplyFromLeft(d_dC).Tr();
  write_file(tmp, tree, "v.bindat");

  // Evaluate bd_plaq (nomenclature as in MILC code)
  tmp = (Pl.val + Pu.val) / 6 / V3;
  tree = tmp.bgf().Tr();
  write_file(tmp, tree, "bd_plaq.bindat");

  // Evaluate st_plaq
  tmp = (Pl.val + Pu.val + Pm.val) / 3 / V3 / T;
  tree = tmp.bgf().Tr();
  write_file(tmp, tree, "st_plaq.bindat");

  // Eval. ss_plaq
  tmp = Ps.val / 3 / V3 / (T-1);
  tree = tmp.bgf().Tr();
  write_file(tmp, tree, "ss_plaq.bindat");
}


struct LocalGaugeFixing1 {
  static const double alpha = .05;
  void operator()(ptGluon_fld& U, const pt::Point& n) const {
    ptSU3 omega;
    for (pt::Direction mu(0); mu.good(); ++mu){
      omega += (U(n, mu) * dag(U (n, mu).bgf()) *
            dag( U(n - mu, mu) ) * U(n - mu, mu).bgf());
    }
    ptSU3 Omega = exp<bgf::AbelianBgf, ORD>( -alpha * omega.reH());
    ptSU3 OmegaDag = exp<bgf::AbelianBgf, ORD>( alpha * omega.reH());
    for (pt::Direction mu(0); mu.good(); ++mu){
      U(n, mu) = Omega * U(n, mu);
      U(n - mu, mu) *= OmegaDag;
    }
  }
};

struct LocalGaugeFixing2 {
  static const double alpha = 0.05;
  void operator()(ptGluon_fld& U, const pt::Point& n) const {
    ptSU3 omega (bgf::unit<bgf::AbelianBgf>());
    for (pt::Direction mu(0); mu.good(); ++mu){
      omega *= (U(n, mu) * dag(U (n, mu).bgf()) *
                dag( U(n + mu, mu) ) * U(n + mu, mu).bgf());
    }
    double atmp = alpha;
    for (int i = 0; i < ORD; ++i){
      W[i] *= atmp;
      atmp *= alpha;
    }
    ptSU3 Winv = inverse(W);
    for (pt::Direction mu(0); mu.good(); ++mu){
      U(n, mu) = W * U(n,mu);
      U(n - mu, mu) *= Winv;
    }  
  }
};

struct LocalGaugeFixing3 {
  static const double alpha = .05;
  void operator()(ptGluon_fld& U, const pt::Point& n) const {
 // exp version
    ptSU3 omega;
    for (pt::Direction mu(0); mu.good(); ++mu)
      omega += U(n, mu) - U(n - mu, mu);

    ptSU3 Omega = exp<bgf::AbelianBgf, ORD>( alpha * omega.reH());
    ptSU3 OmegaDag = exp<bgf::AbelianBgf, ORD>( -alpha * omega.reH());

    for (pt::Direction mu(0); mu.good(); ++mu){
      U(n, mu) = Omega * U(n, mu);
      U(n - mu, mu) *= OmegaDag;
    }
  }
};

struct LocalReUnit {
  void operator()(ptGluon_fld& U, const pt::Point& n) const {
    for (pt::Direction mu(0); mu.good(); ++mu){
      ptsu3 W = get_q(U(n,mu)).reH();
      bgf::AbelianBgf &V = U(n,mu).bgf();
      U(n, mu).ptU() = exp<bgf::AbelianBgf>(W).ptU();
      for (int i = 0; i < ORD; ++i)
        U(n, mu)[i] = V.ApplyFromLeft(U(n,mu)[i]);
    }
  }
};

inline void my_gaugefixing(ptGluon_fld &U){
  static int count = 0;
  int T = act_pars.sz[0] - 1;
  LocalGaugeFixing2 g;
  for (int t = 1; t < T; ++t)
   apply_on_timeslice(U, t, g);
  /*/
  if (++count % 1 == 0){
    LocalReUnit u;
    for (int t = 1; t < T; ++t)
      apply_on_timeslice(U, t, u);
  }
  //*/
}

inline void nu_meas_parallel(const ptGluon_fld &U){
  std::vector<double> nor(ORD*2*3 + 2);
  SU3 C;
  C(0,0) = Cplx(0, 0);
  C(1,1) = Cplx(0, 2./act_pars.sz[1]);
  C(2,2) = Cplx(0, -2./act_pars.sz[1]);


  P_lower El;
  P_upper Eu;
  meas_on_timeslice(U, 0, El);
  meas_on_timeslice(U, act_pars.sz[0] - 2, Eu);
  ptSU3 tmp = El.val + Eu.val;

  for_each(tmp.begin(), tmp.end(), pta::mul(C));
  Cplx tree = tmp.bgf().ApplyFromLeft(C).Tr();
  nor[0] += tree.re;
  nor[1] += tree.im;
  for (int r = 0; r < ORD*2; r += 2){
    for (int i = 0; i < 3; ++i){
      nor[r*3+2*i+2] += tmp[r/2](i,i).re;
      nor[r*3+2*i+3] += tmp[r/2](i,i).im;
    }
  }
  std::ofstream of("nu.txt", std::ios_base::app);
  for (int i = 0; i < ORD*2*3 + 2; ++i)
    of << nor[i] << " ";
  of << std::endl;
  of.close();
}

inline void E_meas_single(ptGluon_fld &Umu) {
  std::vector<double> nor(ORD*2*3 + 2);
  SU3 C;
  const int mu_t = 0;
  C(0,0) = Cplx(0, -1./act_pars.sz[1]);
  C(1,1) = Cplx(0, .5/act_pars.sz[1]);
  C(2,2) = Cplx(0, .5/act_pars.sz[1]);
  ptSU3 tmp;
  for( int x = 0; x < act_pars.sz[1]; ++x)
    for (int y = 0; y < act_pars.sz[2]; ++y)
      for (int z = 0; z < act_pars.sz[3]; ++z){
        //int n = 1 + act_pars.sz[3]*(z + act_pars.sz[2]*(y + act_pars.sz[1]*x) );
        int n = z + act_pars.sz[3]*(y + act_pars.sz[2]*(x + act_pars.sz[1]*1));
        for (int k = 1; k < 4; ++k)
          //   U(n, 0) * U(n + \hat 0, k)
          // * U^\dagger(n + \hat k, 0) * U^\dagger(n, k)
          tmp += 
            Umu.W[Umu.Z->get(n, 1, k, -1, mu_t)][mu_t]*
            dag(Umu.W[Umu.Z->L[n][4]][k])*
            dag(Umu.W[Umu.Z->get(n, -1, mu_t)][mu_t])*
            Umu.W[Umu.Z->get(n, -1, mu_t)][k];
        n = z + act_pars.sz[3]*(y + act_pars.sz[2]*(x + act_pars.sz[1]* (act_pars.sz[0] - 1)));
        for (int k = 1; k < 4; ++k)
          tmp += 
            dag(Umu.W[Umu.Z->L[n][4]][k])*
            dag(Umu.W[Umu.Z->get(n, -1, mu_t)][mu_t])*
            Umu.W[Umu.Z->get(n, -1, mu_t)][k]*
            Umu.W[Umu.Z->get(n, 1, k, -1, mu_t)][mu_t];
          //for_each(tmp.begin(), tmp.end(), pta::mul(C));
          //tree = tmp.bgf().ApplyFromLeft(C).Tr();
          //nor[0] += tree.re;
          //nor[1] += tree.im;
          //for (int r = 0; r < ORD*2; r += 2)
          // for (int i = 0; i < 3; ++i){
          //   nor[r*3+2*i+2] += tmp[r/2](i,i).re;
          //   nor[r*3+2*i+3] += tmp[r/2](i,i).im;
          // }
      }
  for_each(tmp.begin(), tmp.end(), pta::mul(C));
  Cplx tree = tmp.bgf().ApplyFromLeft(C).Tr();
  nor[0] += tree.re;
  nor[1] += tree.im;
  for (int r = 0; r < ORD*2; r += 2){
    for (int i = 0; i < 3; ++i){
      nor[r*3+2*i+2] += tmp[r/2](i,i).re;
      nor[r*3+2*i+3] += tmp[r/2](i,i).im;
    }
  }
  std::ofstream of("E.nor", std::ios_base::app);
  for (int i = 0; i < ORD*2*3 + 2; ++i)
    of << nor[i] << " ";
  of << std::endl;
  of.close();
}

//inline void E_meas(ptGluon_fld &Umu){
//#ifndef __PARALLEL_OMP__
//  E_meas_single(Umu);
//#else
//  E_meas_parallel(Umu);
//#endif
//}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
///
///  check_bgf
///
///  Diagnostics. This function checks if the Background field is
///  modified and throws an exception if so. If you use it, make sure
///  to firt check that do_debug is set to true in newQCDpt.h and try
///  to modify the field, checking if the simulation will indeed
///  terminate.
///
///  \author Dirk Hesse <herr.dirk.hesse@gmail.com>
///  \date Wed Feb  8 15:57:12 2012

inline void check_bgf(ptGluon_fld &Umu) {
  static int j = 0;
  j++;
  for(int i = 0; i < act_pars.iVol; ++i){
    int t = i / (act_pars.sz[1]*act_pars.sz[2]*act_pars.sz[3]); 
    for(int mu = 0; mu < 4; ++mu){
      //if (!(j % 20))
      //    std::cout << "mu = " << mu << ", t = " << t << " :"
      //    << Umu.W[Umu.Z->L[i][4]][mu].bgf() << "  vs.  "
      //    << bgf::get_abelian_bgf(t,mu) << 
      //    "  vs.  " << 
      //      (Umu.W[Umu.Z->L[i][4]][mu].bgf() -bgf::get_abelian_bgf(t, mu))<< std::endl;
      IsZero<true, bgf::AbelianBgf>::check(
       Umu.W[Umu.Z->L[i][4]][mu].bgf() + -bgf::get_abelian_bgf(t, mu));
    }
  // SAFETY CHECK:
  // If you uncomment this, an exception should be thrown.
#ifdef DO_UBER_PARANOID_CHECK
   Umu.W[10][1].bgf() *= 1.00001;
#endif
  }
}


// Beat step di evoluzione quenched
void NsptEvolve(ptGluon_fld& Umu){

#ifdef __TIMING__
  Time.reset();
#endif

  for( int t2 = 0; t2 < nspt_pars.Beat; t2++ ){
    
    // Inizializza a 0 le viariabili di appoggio, le placchette
    // e la norma
    for (int i1=0; i1 < dim; i1++){
      MM[i1].zero();
    }
    for (int i1=0; i1 <= ORD; i1++){
      w1[i1]    = 0.0;
      norm[i1]  = 0.0;
      for (int mu=0; mu < dim; mu++){
	trU[mu][i1]   = 0.0;
      }
#if GAUGE_ACTION != WIL
      w2[i1]    = 0.0;
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

    // Sottrazione dei momenti nulli
#ifdef __TIMING__
    Time.tic_zm();
#endif

    zero_modes_subtraction(Umu);

#ifdef __TIMING__
    Time.toc_zm();
#endif

    // Gauge fixing stocastico
#ifdef __TIMING__
    Time.tic_gf();
#endif

#ifndef __FA_GAUGE_FIXING__
    //stochastic_gauge_fixing(Umu);
    my_gaugefixing(Umu);
#else
    // Fourier acceleration
    FAstochastic_gauge_fixing(Umu);
#endif

#ifdef __TIMING__
    Time.toc_gf();
    Time.reduce();
#endif
    // DH Feb. 7, 2012
#ifdef DO_EXTRA_PARANOID_DEBUG
    // print the norm of U_0 at the t=0 boundary
    //u0_print(Umu);
    // check if the Peter Weisz Background does not get modified.
    check_bgf(Umu);
#endif
  }
  // measure the couplings and various plaquettes
  measure(Umu);

#ifdef __TIMING__
  Time.out();
#endif
}





// Beat step di evoluzione unquenched
void NsptEvolve(ptGluon_fld& Umu, ptSpinColor_fld& Pmu, SpinColor_fld& Xi){

#ifdef __TIMING__
  Time.reset();
#endif

  for( int t2 = 0; t2 < nspt_pars.Beat; t2++ ){
    
    // Inizializza a 0 le viariabili di appoggio, le placchette
    // e la norma
    for (int i1=0; i1 < dim; i1++){
      MM[i1].zero();
    }
    for (int i1=0; i1 <= ORD; i1++){
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

    // Sottrazione dei momenti nulli
#ifdef __TIMING__
    Time.tic_zm();
#endif

    zero_modes_subtraction(Umu);

#ifdef __TIMING__
    Time.toc_zm();
#endif

    // Gauge fixing stocastico
#ifdef __TIMING__
    Time.tic_gf();
#endif

#ifndef __FA_GAUGE_FIXING__
    stochastic_gauge_fixing(Umu);
#else
    // Fourier acceleration
    FAstochastic_gauge_fixing(Umu);
#endif

#ifdef __TIMING__
    Time.toc_zm();
#endif


    // Evoluzione della parte fermionica
#ifdef __TIMING__
    Time.tic_f();
#endif

    fermion_wilson(Umu, Pmu, Xi);

#ifdef __TIMING__
    Time.toc_f();
    Time.reduce();
#endif

  }

#ifdef __TIMING__
  Time.out();
#endif

}



// Confronta 2 placchette. Se sono uguali ritorna 0, altrimenti
// la differenza tra le due
int plaquette_check(Cplx* w1, Cplx* w2){
  // return ( memcmp(w1, w2, (ORD+1)*sizeof(Cplx)) ); 
  for(int i1 = 0; i1 <= ORD;i1++){
    if( (w1[i1]-w2[i1]).mod() > 1e-7 ) return 1;
  }
  return 0;
}





#ifdef __WILLOOP_2x2__

void WL2x2(ptGluon_fld& Umu) {
  
#ifndef __PARALLEL_OMP__
  ptSU3 tmp, tmp1;
#else
  ptSU3 tmp[NTHR], tmp1[NTHR];
#endif
  
  int a,b,d,e;
  

#ifndef __PARALLEL_OMP__
  for(int site_x = 0; site_x < act_pars.iVol; site_x++) {

    for(int mu = 3; mu > 0; --mu){
      a = Umu.Z->L[site_x][5+mu];
      b = Umu.Z->L[a][5+mu];
      tmp = Umu.W[site_x][mu] * Umu.W[a][mu];

      for( int nu = 0; nu < mu; ++nu){
	d = Umu.Z->L[site_x][5+nu];
	e = Umu.Z->L[d][5+nu];
	tmp1 += tmp * ( Umu.W[b][nu] * Umu.W[ Umu.Z->L[b][5+nu] ][nu] * 
			dag( Umu.W[site_x][nu] * Umu.W[d][nu] * Umu.W[e][mu] * 
			     Umu.W[ Umu.Z->L[e][5+mu] ][mu] ));

      } // nu
    } // mu
  } // site

  for( int i1 = 1; i1 < ORD; i1+=2){
    loop2x2 << (tmp1[i1].Tr()).re*act_pars.rVol*D18 << "\t";
  }
  loop2x2 << std::endl;

#else
#pragma omp parallel private(tid,a,b,d,e)  num_threads(NTHR) 
  {
    int chunk = act_pars.iVol/NTHR;
    tid = omp_get_thread_num();
    for(int site_x = chunk*tid; site_x < chunk*(tid+1);site_x++){
      
      for(int mu = 3; mu > 0; --mu){
	a = Umu.Z->L[site_x][5+mu];
	b = Umu.Z->L[a][5+mu];
	tmp[tid] = Umu.W[site_x][mu] * Umu.W[a][mu];
	
	for( int nu = 0; nu < mu; ++nu){
	  d = Umu.Z->L[site_x][5+nu];
	  e = Umu.Z->L[d][5+nu];
	  tmp1[tid] += tmp[tid] * ( Umu.W[b][nu] * Umu.W[ Umu.Z->L[b][5+nu] ][nu] * 
				    dag( Umu.W[site_x][nu] * Umu.W[d][nu] * Umu.W[e][mu] * 
					 Umu.W[ Umu.Z->L[e][5+mu] ][mu] ));
	} // nu
      } // mu
      
    } // end site_x
    
  } // end parallel
  

  // Reduce
  for(int tid = 1; tid < NTHR; tid++){
    tmp1[0] += tmp1[tid];
  }

  // std::cout << std::endl;
  for( int i1 = 1; i1 <= ORD; i1+=2){
    // std::cout << (tmp1[0].ptU[i1].Tr()).re*act_pars.rVol*D18 << "\t";
    loop2x2 << (tmp1[0][i1].Tr()).re*act_pars.rVol*D18 << "\t";
  }

  // std::cout << std::endl;
  loop2x2 << std::endl;

#endif

}

#endif


#ifdef __MANY_LOOP__

Cplx l21[ORD+1];
Cplx l22[ORD+1];
Cplx l42[ORD+1];
Cplx l44[ORD+1];


void ComputeLoops(ptGluon_fld &Umu){

  ptSU3 br, bl, mr, ml, tr, tl;
  ptSU3 ld, lu, cd, cu, rd, ru;


  ptSU3 w22,w21,w44,w42;
#ifdef __PARALLEL_OMP__
  ptSU3 t_w22[NTHR],t_w21[NTHR],t_w44[NTHR],t_w42[NTHR];
#endif

  ptSU3 tmp;

  int a,b,c,d,e,f,g,h,j,k;
  int i,l,m,n,o,p,q,r,s,t;

  w22.zero();
  w21.zero();
  w42.zero();
  w44.zero();

#ifdef __PARALLEL_OMP__
  for( int thr = 0; thr < NTHR; thr++)
    {
      t_w22[thr].zero();
      t_w21[thr].zero();
      t_w42[thr].zero();
      t_w44[thr].zero();    
    }
#endif

  for( int mu = 0; mu < dim-1; mu++)
    {
      for( int nu = mu+1; nu < dim; nu++)
	{

#ifdef __PARALLEL_OMP__
#pragma omp parallel for private(a,b,c,d,e,f,g,h,j,k, i,l,m,n,o,p,q,r,s,t,br, bl, mr, ml, tr, tl, ld, lu, cd, cu, rd, ru) num_threads(NTHR)
#endif
	  for( int x = 0; x < Umu.Z->Size; x++)
	    {
#ifdef __PARALLEL_OMP__
	      int thr = omp_get_thread_num();
#endif
	      a = Umu.Z->get(x, 1, mu);
	      b = Umu.Z->get(a, 1, mu);
	      c = Umu.Z->get(b, 1, nu);
	      e = Umu.Z->get(x, 1, nu);
	      d = Umu.Z->get(e, 1, mu);

  	      f = Umu.Z->get(c, 1, nu);
	      g = Umu.Z->get(d, 1, nu);
	      h = Umu.Z->get(e, 1, nu);


	      bl = Umu.W[x][mu]*Umu.W[a][mu];
	      ml = Umu.W[h][mu]*Umu.W[g][mu];
	      cd = Umu.W[b][nu]*Umu.W[c][nu];
	      ld = Umu.W[x][nu]*Umu.W[e][nu];

#ifdef  __PARALLEL_OMP__
	      t_w22[thr] += bl*cd*dag(ld*ml);
	      
	      t_w21[thr] += ( bl*Umu.W[b][nu]*
			      dag(Umu.W[x][nu]*Umu.W[e][mu]*Umu.W[d][mu]) +
			      Umu.W[x][mu]*Umu.W[a][nu]*Umu.W[d][nu]*
			      dag( ld*Umu.W[h][mu]) );	      
#else
	      w22 += bl*cd*dag(ld*ml);

	      w21 += ( bl*Umu.W[b][nu]*
		       dag(Umu.W[x][nu]*Umu.W[e][mu]*Umu.W[d][mu]) +
		       Umu.W[x][mu]*Umu.W[a][nu]*Umu.W[d][nu]*
		       dag( ld*Umu.W[h][mu]) );	      
#endif

              i = Umu.Z->get(h, 1, nu);
              j = Umu.Z->get(f, 1, nu);
              k = Umu.Z->get(j, 1, nu);
              n = Umu.Z->get(i, 1, nu);
              l = Umu.Z->get(n, 1, mu);
              m = Umu.Z->get(b, 1, mu);
              o = Umu.Z->get(m, 1, mu);
              p = Umu.Z->get(o, 1, nu);
              q = Umu.Z->get(p, 1, nu);
              r = Umu.Z->get(f, 1, mu);
              s = Umu.Z->get(q, 1, nu);
              t = Umu.Z->get(k, 1, mu);

              br = Umu.W[b][mu]*Umu.W[m][mu];
              mr = Umu.W[f][mu]*Umu.W[r][mu];
              tr = Umu.W[k][mu]*Umu.W[t][mu];
              tl = Umu.W[n][mu]*Umu.W[l][mu];
              lu = Umu.W[h][nu]*Umu.W[i][nu];             
              cu = Umu.W[f][nu]*Umu.W[j][nu];
              ru = Umu.W[q][nu]*Umu.W[s][nu];
              rd = Umu.W[o][nu]*Umu.W[p][nu];
              

#ifdef  __PARALLEL_OMP__
              t_w42[thr] += bl*(br*rd*dag(ml*mr) 
				+ cd*cu*dag(lu*tl))*dag(ld);
	      
              t_w44[thr] += bl*br*rd*ru*dag(ld*lu*tl*tr);	      
#else
              w42 += bl*(br*rd*dag(ml*mr) 
                         + cd*cu*dag(lu*tl))*dag(ld);

              w44 += bl*br*rd*ru*dag(ld*lu*tl*tr);
#endif
	    } // i
	} // nu
    } // mu



#ifdef  __PARALLEL_OMP__
  for( int thr = 0; thr < NTHR; thr++)
    {
      w21+=t_w21[thr];
      w22+=t_w22[thr];
      w42+=t_w42[thr];
      w44+=t_w44[thr];
    }
#endif


  w21.Tr(l21);
  w22.Tr(l22);
  w42.Tr(l42);
  w44.Tr(l44);

  const double loopNorm = 1./(3*2*3*(Umu.Z->Size));
  // loop1x2 << "Plaquette 1x2" << endl;
  for(int i1 = 0; i1 <= ORD; i1+=2){  
    loop2x2 << l21[i1].re*.5*loopNorm << "\t";
  }					
  
  // loop2x2 << "Plaquette 2x2" << endl;
  for(int i1 = 0; i1 <= ORD; i1+=2){  
    loop2x2 << l22[i1].re*loopNorm << "\t";
  }					
  
  // loop2x2 << "Plaquette 2x4" << endl;
  for(int i1 = 0; i1 <= ORD; i1+=2){  
    loop2x2 << l42[i1].re*.5*loopNorm << "\t";
  }					
  
  // cout << "Plaquette 4x4" << endl;
  for(int i1 = 0; i1 <= ORD; i1+=2){  
    loop2x2 << l44[i1].re*loopNorm << "\t";
  }					
  loop2x2 << std::endl;

}

#endif


// Misura a configurazione congelata la placchetta
void plaquette_measure(ptGluon_fld& Umu, act_params_t& act_pars){
#ifndef __PARALLEL_OMP__
  for(int i1 = 0; i1 < ORD+1; i1++){
    w[i1] = 0.0;
  }					
  
  Cplx* app = new Cplx[ORD+1];
  
  for(int i = 0; i < act_pars.iVol;i++){
    W1.zero();				
    link_c = Umu.Z->L[i][dim];
    for(int mu = 0; mu < dim-1; mu++){
      W.zero();				
      for(int nu = mu+1; nu < dim; nu++){
	if(nu != mu ){
	  W += Umu.staple(i, mu, nu);
	}
      }
      W1 = Umu.W[link_c][mu]*W;
      W1.Tr(app);
      for(int i1 = 0; i1 <= ORD; i1++){
	w[i1] += app[i1];
      }
    }	
  }
     
  for(int i1 = 0; i1 <= ORD; i1++){
    w[i1] /= (act_pars.iVol*36);
  }					

  W.zero();
  W1.zero();
  delete [] app;

#else

  for(int nt = 0; nt < NTHR; nt++){
    for(int i1 = 0; i1 < ORD+1; i1++){
      ww[nt][i1] = 0.0;
    }				
  }
  //  Cplx* app = new Cplx[ORD+1];
  
#pragma omp parallel private(tid,link_c) num_threads(NTHR)
  {
    int chunk = act_pars.iVol/NTHR;

    tid = omp_get_thread_num();    
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
	// 	Ww1[tid] = Umu.W[link_c][mu]*Ww[tid];
 	Ww1[tid] = Umu.W[site_x][mu]*Ww[tid];
 	Ww1[tid].Tr(ww1[tid]);
	
 	for(int i1 = 0; i1 <= ORD; i1++){
 	  ww[tid][i1] += ww1[tid][i1];
 	}
	
      }// altro link, nello stesso sito	
    } // fine siti
    
    Ww[tid].zero();
    Ww1[tid].zero();
  } // end parallel
#pragma omp barrier

  for(int i1 = 0; i1 <= ORD; i1++){
    w[i1] = 0;
    for(int nt = 0; nt < NTHR; nt++){
      //      std::cout << ww[nt][i1].re << "\t";
      w[i1] += ww[nt][i1];
    }
    //    std::cout << std::endl;
    w[i1] /= (act_pars.iVol*72);		
  }					
  //  std::cout << std::endl;
#endif

}


// // Misura a configurazione congelata la placchetta
// void plaquette_measure(ptGluon_fld& Umu, act_params_t& act_pars){
// #ifndef __PARALLEL_OMP__
//   for(int i1 = 0; i1 < ORD+1; i1++){
//     w[i1] = 0.0;
//   }					
  
//   Cplx* app = new Cplx[ORD+1];
  
//   for(int i = 0; i < act_pars.iVol;i++){
//     W1.zero();				
//     link_c = Umu.Z->L[i][dim];
//     for(int mu = 0; mu < dim; mu++){
//       W.zero();				
//       for(int nu = 0; nu < dim; nu++){
// 	if(nu != mu ){
// 	  W += Umu.staple(i, mu, nu);
// 	}
//       }
//       W1 = Umu.W[link_c][mu]*W;
//       W1.Tr(app);
//       for(int i1 = 0; i1 <= ORD; i1++){
// 	w[i1] += app[i1];
//       }
//     }	
//   }
     
//   for(int i1 = 0; i1 <= ORD; i1++){
//     w[i1] /= (act_pars.iVol*72);
//   }					

//   W.zero();
//   W1.zero();
//   delete [] app;

// #else

//   for(int nt = 0; nt < NTHR; nt++){
//     for(int i1 = 0; i1 < ORD+1; i1++){
//       ww[nt][i1] = 0.0;
//     }				
//   }
//   //  Cplx* app = new Cplx[ORD+1];
  
// #pragma omp parallel private(tid,link_c) num_threads(NTHR)
//   {
//     tid = omp_get_thread_num();
    
//     for(int site_x = chunk*tid; site_x < chunk*(tid+1);site_x++){
//       Ww1[tid].zero();
//       //      link_c = Umu.Z->L[site_x][dim];
//       for(int mu = 0; mu < dim; mu++){
//  	Ww[tid].zero();
//  	for(int nu = 0; nu < dim; nu++){
//  	  if(nu != mu ){
//  	    Ww[tid] += Umu.staple(site_x, mu, nu);
//  	  }
//  	}
// 	// 	Ww1[tid] = Umu.W[link_c][mu]*Ww[tid];
//  	Ww1[tid] = Umu.W[site_x][mu]*Ww[tid];
//  	Ww1[tid].Tr(ww1[tid]);
	
//  	for(int i1 = 0; i1 <= ORD; i1++){
//  	  ww[tid][i1] += ww1[tid][i1];
//  	}
	
//       }// altro link, nello stesso sito	
//     } // fine siti
    
//     Ww[tid].zero();
//     Ww1[tid].zero();
//   } // end parallel
// #pragma omp barrier

//   for(int i1 = 0; i1 <= ORD; i1++){
//     w[i1] = 0;
//     for(int nt = 0; nt < NTHR; nt++){
//       //      std::cout << ww[nt][i1].re << "\t";
//       w[i1] += ww[nt][i1];
//     }
//     //    std::cout << std::endl;
//     w[i1] /= (act_pars.iVol*72);		
//   }					
//   //  std::cout << std::endl;
//  #endif

// }