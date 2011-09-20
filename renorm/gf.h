#include "../QCDenvNODEpt.h"
#include "renorm_pars.h"

#define GF_PREC 0.001

/********************** MACRO **************************/

#define PLAQUETTE_MEAS(a)				\
  {							\
    for(int i1 = 0; i1 < PTORD+1; i1++){		\
      a[i1]  = 0.0;					\
    }							\
    for(int i = 0; i < iVol; i++){			\
      for(int mu = 0; mu < dim; mu++){			\
	F.zero();					\
	for(int nu = 0; nu < dim; nu++){		\
	  if(nu != mu){					\
	    F += Umu.staple(i, mu, nu);			\
	  }						\
	}						\
	FF += U[get(&Umu, i, mu)]*F;			\
	FF.Tr(a);					\
      }							\
    }							\
    F.zero();						\
    FF.zero();						\
    for(int i1 = 0; i1 < PTORD+1; i1++){		\
      a[i1] /= (iVol*72.0);				\
    }							\
  }


#define PLAQUETTE_CHECK(a,b){						\
    if(memcmp(a,b,sizeof(Cplx)*(PTORD+1))){				\
      printf("La misura della placchetta non corrisponde\n");		\
      exit(1);								\
    }									\
    memset(a,0,sizeof(Cplx)*(PTORD+1));					\
    memset(b,0,sizeof(Cplx)*(PTORD+1));					\
  }

/********************** FINE MACRO **************************/


void dirac(int, int, ptSpinColor_fld*, ptGluon_fld&);
void M0invT(ptSpinColor_fld*);
void diracT(int ptord, ptGluon_fld&, ptSpinColor_fld*);


using namespace std;

int carica_campo(ptGluon_fld &Umu){

  Cplx *plaq, *plaq2;
  ptSU3 *U;
  ptSU3 F,FF;
  int iVol = Umu.Z->Size;

  plaq  = new Cplx[PTORD+1];
  plaq2 = new Cplx[PTORD+1];

  Umu.load(conf, plaq);
  U = Umu.handle();

  PLAQUETTE_MEAS(plaq2);
  PLAQUETTE_CHECK(plaq, plaq2);

  delete[] plaq;
  delete[] plaq2;
  return 0;
}




int gauge_fixing(latt &LL, ptGluon_fld &Umu, double &alpha, int &ptord, double &old, int &cnt){

#ifdef ABC
  LL.p_pbc();
#endif
  ptGluon_fld copia(&LL);
  ptSU3_fld W(&LL);

  ptSU3 *U,zero;
  U = Umu.handle();

  Cplx accum[ptord+1];
  int coord[4];

  accum[ptord-1].re = 1;
  old = 10000*LL.Size;
  copia = Umu;
  
  
  fftw_plan plan[2];
  cout << "\nAlloco i plan per le fft..." ;
  plan[0] = fftw_plan_many_dft(dim, LL.Sz, NC*NC*PTORD,
			       (fftw_complex *) W.W, NULL, NC*NC*PTORD+1, 1,
			       (fftw_complex *) W.W, NULL, NC*NC*PTORD+1, 1,
			       FFTW_FORWARD, FFTW_MEASURE);
  plan[1] = fftw_plan_many_dft(dim, LL.Sz, NC*NC*PTORD,
			       (fftw_complex *) W.W, NULL, NC*NC*PTORD+1, 1,
			       (fftw_complex *) W.W, NULL, NC*NC*PTORD+1, 1,
			       FFTW_BACKWARD, FFTW_MEASURE);
  cout << "fatto\n" 
       << "Fisso la gauge...";

  while(fabs(accum[ptord-1].re) > GF_PREC){
    cnt++;
    
    memset(accum,0,2*(ptord+1)*sizeof(double));
    
    for(int x = 0; x < LL.Size; x++){
      W.W[x] = zero;
      for(int mu = 0; mu < dim; mu++){
	W.W[x] += log(Umu.W[x].U[mu]) - log(U[get(&Umu,x,-1,mu,mu)]);
      }
      // norma campo W(x)
      for(int o = 0; o < ptord; o++){
	accum[o] += (W.W[x].ptU[o]*dag(W.W[x]).ptU[o]).Tr();
      }
    }
    cout << "Norma: " << accum[ptord-1].re << endl;

    
    //    accelerazione di Fourier
    fftw_execute(plan[0]);
    for(int x = 1; x < LL.Size; x++){
      LL.get(x,coord);
      W.W[x] /= -((LL.p2hat[0][coord[0]]+
    		  LL.p2hat[1][coord[1]]+
    		  LL.p2hat[2][coord[2]]+
    		  LL.p2hat[3][coord[3]])
    		 *LL.Size);
    }
    fftw_execute(plan[1]);
    
    
    for(int x = 0; x < LL.Size; x++){
      W.W[x] = exp(alpha*W.W[x]);
    }
    
    for(int x = 0; x < LL.Size; x++){
      for(int mu = 0; mu < dim; mu++){
	U[get(&Umu,x,mu)] = W.W[x]*U[get(&Umu,x,mu)]*dag(W.get(x,1,mu));
      }
    }

    
    if(old < accum[ptord-1].re){
/*       old = 1./0; */
      Umu = copia;
      return 1;
    }
    old = accum[ptord-1].re;
  }
  cout << "fatto!" << endl << endl;

  fftw_destroy_plan(plan[0]);
  fftw_destroy_plan(plan[1]);

#ifdef ABC
  LL.p_abc();
#endif

  return 0;
}



void Minv(ptSpinColor_fld *Pmu, int ptmax, ptGluon_fld &Umu, int tra = 0){
  // parto da una delta in p, applico l'ordine 0 in p.

  if(tra == 0){
    Pmu->M0inv();
  }
  else{
    M0invT(Pmu);
  }
  // passo in q, e copio nel campo spinore perturbativo
  Pmu->fft(1);
  
  for(int i = 0; i < Pmu->Z->Size; i++){
    for(int mu = 0; mu < dim; mu++){
      Pmu->psi[i].psi[mu].ptCV[0] = Pmu->scfld->psi[i].psi[mu]; // resta in Q
    }
  }
  
  for(int ord = 1; ord <= ptmax; ord++){    

    // questa chiamata sfrutta direttamente gli ordini perturbativi in memoria
    // che devono essere nello spazio delle posizione!
    if(tra == 0){
      Pmu->dirac(ord, Umu);
      
      // lo passo in p e applico M0^(-1)
      Pmu->fft(0);
      Pmu->M0inv();
    }
    else{
      diracT(ord, Umu, Pmu);
      
      // lo passo in p e applico M0^(-1)
      Pmu->fft(0);
      M0invT(Pmu);
    }
    // prima di mettere nel perturbativo torno alle posizioni
    Pmu->fft(1);
    
    for(int i = 0; i < Pmu->Z->Size; i++){
      for(int mu = 0; mu < dim; mu++){
	Pmu->psi[i].psi[mu].ptCV[ord] = -(Pmu->scfld->psi[i].psi[mu]);
      }
    }
  }

  // alla fine devo rimettere in p (passo attraverso scfld)
  for(int ord = 0; ord <= ptmax; ord++){    
    for(int i = 0; i < Pmu->Z->Size; i++){
      for(int mu = 0; mu < dim; mu++){
	Pmu->scfld->psi[i].psi[mu] = Pmu->psi[i].psi[mu].ptCV[ord];
      }
    }
    Pmu->fft(0);
    for(int i = 0; i < Pmu->Z->Size; i++){
      for(int mu = 0; mu < dim; mu++){
	Pmu->psi[i].psi[mu].ptCV[ord] = Pmu->scfld->psi[i].psi[mu];
      }
    }
  }

}





void M0invT(ptSpinColor_fld *Pmu){
  double gmT[4] = {1.0, -1.0, 1.0, -1.0};
  
  double diag, den;
  Cplx im_(0, -1);
  int i = 0, pp[4];
  
#ifdef ABC    
  Pmu->scfld->psi[0] *= Cplx(0.0,0.0);
    for (pp[0] = 0; pp[0] < Pmu->Z->Sz[0]; pp[0]++) {
      for (pp[1] = 0; pp[1] < Pmu->Z->Sz[1]; pp[1]++) {
	for (pp[2] = 0; pp[2] < Pmu->Z->Sz[2]; pp[2]++) {		
	  for (pp[3] = 0; pp[3] < Pmu->Z->Sz[3]; pp[3]++) {
	    i = Pmu->Z->get(pp);
	    diag = MBARE + 0.5*(Pmu->Z->p2hat[0][pp[0]] + Pmu->Z->p2hat[1][pp[1]] + 
				Pmu->Z->p2hat[2][pp[2]] + Pmu->Z->p2hat[3][pp[3]]);
	    den = diag*diag + (Pmu->Z->p2bar[0][pp[0]] + Pmu->Z->p2bar[1][pp[1]] + 
			       Pmu->Z->p2bar[2][pp[2]] + Pmu->Z->p2bar[3][pp[3]]);
	    
	    Pmu->scfld->psi[i] = ( ((Pmu->scfld->psi[i]).gmleft(0)*gmT[0] * (Pmu->Z->pbar[0][pp[0]]) +
				    (Pmu->scfld->psi[i]).gmleft(1)*gmT[1] * (Pmu->Z->pbar[1][pp[1]]) +
				    (Pmu->scfld->psi[i]).gmleft(2)*gmT[2] * (Pmu->Z->pbar[2][pp[2]]) +
				    (Pmu->scfld->psi[i]).gmleft(3)*gmT[3] * (Pmu->Z->pbar[3][pp[3]]))
				   *im_+ (Pmu->scfld->psi[i])*diag) / den;
	  }
	}
      }
    }
#elif defined PBC
    Pmu->scfld->psi[0] *= Cplx(0.0,0.0);
    for (pp[0] = 0; pp[0] < Pmu->Z->Sz[0]; pp[0]++) {
      for (pp[1] = 0; pp[1] < Pmu->Z->Sz[1]; pp[1]++) {
	for (pp[2] = 0; pp[2] < Pmu->Z->Sz[2]; pp[2]++) {		
	  for (pp[3] = 0; pp[3] < Pmu->Z->Sz[3]; pp[3]++) {
	    if( (pp[0]+pp[1]+pp[2]+pp[3])!= 0){
	      i = Pmu->Z->get(pp);
	      diag = MBARE + 0.5*(Pmu->Z->p2hat[0][pp[0]] + Pmu->Z->p2hat[1][pp[1]] + 
				  Pmu->Z->p2hat[2][pp[2]] + Pmu->Z->p2hat[3][pp[3]]);
	      den = diag*diag + (Pmu->Z->p2bar[0][pp[0]] + Pmu->Z->p2bar[1][pp[1]] + 
				 Pmu->Z->p2bar[2][pp[2]] + Pmu->Z->p2bar[3][pp[3]]);
	      
	      Pmu->scfld->psi[i] = ( ((Pmu->scfld->psi[i]).gmleft(0)*gmT[0] * (Pmu->Z->pbar[0][pp[0]]) +
				      (Pmu->scfld->psi[i]).gmleft(1)*gmT[1] * (Pmu->Z->pbar[1][pp[1]]) +
				      (Pmu->scfld->psi[i]).gmleft(2)*gmT[2] * (Pmu->Z->pbar[2][pp[2]]) +
				      (Pmu->scfld->psi[i]).gmleft(3)*gmT[3] * (Pmu->Z->pbar[3][pp[3]]))
				     *im_+ (Pmu->scfld->psi[i])*diag) / den;
	    }
	  }
	}
      }
    }
#endif
  }


void diracT(int ptord, ptGluon_fld &Umu, ptSpinColor_fld *Pmu)  {
  SpinColor S1, S2, S;
  double gmT[4] = {1.0, -1.0, 1.0, -1.0};
  
  for (int i = 0; i < Pmu->Z->Size; i++) {
    S *= 0;	
    
    for (int ord = ptord; ord > 0; ord--) {
      
      for (int nu = 0; nu < 4; nu++) {
	S.psi[nu] += ( Pmu->get(i).psi[nu].ptCV[ptord-ord] * mcpt[ord] );
      }
      
      for (int mu = 0; mu < 4; mu++) {
	
	for (int nu = 0; nu < 4; nu++) {
	  S1.psi[nu] = ~(Umu.get(i, mu).ptU[ord-1])*
	    Pmu->near(i, +1, mu).psi[nu].ptCV[ptord-ord];
	  S2.psi[nu] = tra(Umu.near(i, -1, mu).U[mu].ptU[ord-1])*
	    Pmu->near(i, -1, mu).psi[nu].ptCV[ptord-ord];
	}
	
	S1 -= S1.gmleft(mu)*gmT[mu];
	S2 += S2.gmleft(mu)*gmT[mu];	  
	S -= (S1 + S2)*.5;
      }
    }
    
    Pmu->scfld->psi[i] = S;
    
  }
}
