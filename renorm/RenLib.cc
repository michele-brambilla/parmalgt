#include "RenLib.h"
#include <GammaStuff.h>
#include <MassC.h>


extern renorm_params_t params;

using namespace std;

int carica_campo(ptGluon_fld &Umu, Cplx* plaq){

  if( Umu.load(params.conf, plaq) ){
    //  if(!Umu.load_ape(conf)){
     cout << "Configurazione non trovata." << endl;
     return 1;
   }

  return 0;
}




int gauge_fixing(ptGluon_fld &Umu, int &cnt){

  ptSU3_fld W(Umu.Z);
  SU3_fld W1(Umu.Z);

  ptSU3 *U,zero;
  U = Umu.handle();

  Cplx accum[allocORD+1];
  int coord[4];

  accum[params.ptord].re = 1;
  double old = 2.0e+23;
  
#ifdef __RENORM_OMP__
  if( !fftw_init_threads() )
    {
      cout << "FFTW multithread not available, single thread execution." 
	   << endl;
    }
  else
    {
      fftw_plan_with_nthreads(NTHR);
    }
#endif  
  
  fftw_plan plan[2];
  cout << "\nAlloco i plan per le fft..." ;
   

  plan[0] = fftw_plan_many_dft(dim, params.sz, NC*NC*PTORD+1,
			       (fftw_complex *) W.W, NULL, NC*NC*PTORD+1, 1,
			       (fftw_complex *) W.W, NULL, NC*NC*PTORD+1, 1,
			       FFTW_FORWARD, FFTW_ESTIMATE);
  plan[1] = fftw_plan_many_dft(dim, params.sz, NC*NC*PTORD+1,
			       (fftw_complex *) W.W, NULL, NC*NC*PTORD+1, 1,
			       (fftw_complex *) W.W, NULL, NC*NC*PTORD+1, 1,
			       FFTW_BACKWARD, FFTW_ESTIMATE);

  cout << "fatto\n" 
       << "Fisso la gauge...";

  ofstream logfile;
  char *logname;
  time_t rawtime;
  struct tm * timeinfo;
  
  time (&rawtime );
  timeinfo = localtime ( &rawtime );

  logname = new char[100];

  unsigned int i = 0;
  while(i < strlen(params.conf))
    {
#if GAUGE_ACTION == WIL
      if(strncmp(params.conf+i,"wil",3) == 0) break;   
#elif GAUGE_ACTION == TLSYM
      if(strncmp(params.conf+i,"tls",3) == 0) break;
#elif GAUGE_ACTION == IWA
      if(strncmp(params.conf+i,"iwa",3) == 0) break;
#elif GAUGE_ACTION == DBW2
      if(strncmp(params.conf+i,"dbw",3) == 0) break;
#endif
      i++;
    }
  strcpy(logname,"gf_");
  strcat(logname,params.conf+i);

  cout << "See log: " << logname << endl;
  logfile.open(logname);

  logfile << "configurazione:\t" << params.conf    << endl;
  logfile << "precisione:\t"     << GF_PREC        << endl;
  logfile << "ordine:\t"         << params.ptord-1 << endl;

  logfile << "Gauge fixing started at " << ctime(&rawtime) << endl << endl;
  
  while(old > GF_PREC){
    cnt++;
    
    for(int o = 0; o < params.ptord; o++){
      accum[o] = 0;
    }    
    
    for(int x = 0; x < params.iVol; x++){
      W.W[Umu.Z->get(x)] = zero;
      for(int mu = 0; mu < dim; mu++){
	W.W[Umu.Z->get(x)] += log(U[get(&Umu,x,mu)]) - log(U[get(&Umu,x,-1,mu,mu)]);
	//	W.W[LL.get(x)] += U[get(&Umu,x,mu)] - U[get(&Umu,x,-1,mu,mu)];
      }

      W.W[Umu.Z->get(x)] -= dag(W.W[Umu.Z->get(x)]);
      // norma campo W(x)      
      for(int o = 0; o < params.ptord; o++)
	accum[o] += (W.W[x].ptU[o]*dag(W.W[x].ptU[o])).Tr();
    }
    for(int o = 0; o < params.ptord; o++)
      logfile << accum[o].mod() << endl;
    logfile << endl;
    logfile.flush();
    
    fftw_execute(plan[0]);
    
    W.W[0] *= 0.0;
       
#ifdef __RENORM_OMP__
#pragma omp for schedule(static)
#endif
    for(int x = 1; x < Umu.Z->Size; x++){
      Umu.Z->get(x,coord);
      W.W[Umu.Z->get(x)] *= -(params.alpha/(Umu.Z->p2hat[0][coord[0]]+
					    Umu.Z->p2hat[1][coord[1]]+
					    Umu.Z->p2hat[2][coord[2]]+
					    Umu.Z->p2hat[3][coord[3]]));
    }
  
  
    fftw_execute(plan[1]);
  
  
#ifdef __RENORM_OMP__
#pragma omp for schedule(static)
#endif
    for(int x = 0; x < Umu.Z->Size; x++){
      W.W[x] /= (double)Umu.Z->Size;
      W.W[x].reH();
      W.W[x] = exp(W.W[x]);
    }
  
  
#ifdef __RENORM_OMP__
#pragma omp for schedule(static)
#endif
    for(int x = 0; x < Umu.Z->Size; x++){
      for(int mu = 0; mu < dim; mu++){
	U[get(&Umu,x,mu)] = W.get(x)*Umu.get(x,mu)*dag(W.get(x,1,mu));
      }
    }
    
    if(old < accum[params.ptord-1].mod()){
      //if(old < accum[ptord-1].re){
      /*       Umu = copia; */
      return 1;
    }
    old = accum[params.ptord-1].mod();
    //old = accum[ptord-1].re;
  }
  cout << "fatto!" << endl << endl;
  

  fftw_destroy_plan(plan[0]);
  fftw_destroy_plan(plan[1]);

#ifdef _FFTW_THR
  fftw_cleanup_threads();
#endif


  logfile << endl << "numero iterazioni = " << cnt << endl;

  time (&rawtime );
  logfile << endl << "Gauge fixing ended at " << ctime(&rawtime) << endl << endl;

  logfile.close();
  delete [] logname;

  return 0;
}



void fft(ptSpinColor_fld** Psi, int ind, int ord) {

#ifndef __RENORM_OMP__
  for( int i = 0 ; i < Psi[0]->Z->Size; i++ ) {
    for (int mu = 0; mu < dim; mu++ ) {
      Psi[0]->scfld->psi[i].psi[mu] = Psi[0]->psi[i].psi[mu].ptCV[ord];
      Psi[1]->scfld->psi[i].psi[mu] = Psi[1]->psi[i].psi[mu].ptCV[ord];
      Psi[2]->scfld->psi[i].psi[mu] = Psi[2]->psi[i].psi[mu].ptCV[ord];
      Psi[3]->scfld->psi[i].psi[mu] = Psi[3]->psi[i].psi[mu].ptCV[ord];
      Psi[4]->scfld->psi[i].psi[mu] = Psi[4]->psi[i].psi[mu].ptCV[ord];
      Psi[5]->scfld->psi[i].psi[mu] = Psi[5]->psi[i].psi[mu].ptCV[ord];
      Psi[6]->scfld->psi[i].psi[mu] = Psi[6]->psi[i].psi[mu].ptCV[ord];
      Psi[7]->scfld->psi[i].psi[mu] = Psi[7]->psi[i].psi[mu].ptCV[ord];
    }
  }
  
  Psi[0]->fftT(ind);
  Psi[1]->fftT(ind);
  Psi[2]->fftT(ind);
  Psi[3]->fftT(ind);
  Psi[4]->fft(ind);
  Psi[5]->fft(ind);
  Psi[6]->fft(ind);
  Psi[7]->fft(ind);

  for( int i = 0 ; i < Psi[0]->Z->Size; i++ ) {
    for (int mu = 0; mu < dim; mu++ ) {
      Psi[0]->psi[i].psi[mu].ptCV[ord] = Psi[0]->scfld->psi[i].psi[mu];
      Psi[1]->psi[i].psi[mu].ptCV[ord] = Psi[1]->scfld->psi[i].psi[mu];
      Psi[2]->psi[i].psi[mu].ptCV[ord] = Psi[2]->scfld->psi[i].psi[mu];
      Psi[3]->psi[i].psi[mu].ptCV[ord] = Psi[3]->scfld->psi[i].psi[mu];
      Psi[4]->psi[i].psi[mu].ptCV[ord] = Psi[4]->scfld->psi[i].psi[mu];
      Psi[5]->psi[i].psi[mu].ptCV[ord] = Psi[5]->scfld->psi[i].psi[mu];
      Psi[6]->psi[i].psi[mu].ptCV[ord] = Psi[6]->scfld->psi[i].psi[mu];
      Psi[7]->psi[i].psi[mu].ptCV[ord] = Psi[7]->scfld->psi[i].psi[mu];
    }
  }
#else
#pragma omp parallel num_threads(NTHR)  shared(Psi)
  {
    int tid = omp_get_thread_num();
    for( int i = 0 ; i < Psi[0]->Z->Size; i++ ) {
      for (int mu = 0; mu < dim; mu++ ) {
	Psi[tid    ]->scfld->psi[i].psi[mu] = Psi[tid    ]->psi[i].psi[mu].ptCV[ord];
	Psi[tid+dim]->scfld->psi[i].psi[mu] = Psi[tid+dim]->psi[i].psi[mu].ptCV[ord];
      }
    }
  
    Psi[tid    ]->fftT(ind);
    Psi[tid+dim]->fft(ind);

    for( int i = 0 ; i < Psi[0]->Z->Size; i++ ) {
      for (int mu = 0; mu < dim; mu++ ) {
	Psi[tid    ]->psi[i].psi[mu].ptCV[ord] = Psi[tid    ]->scfld->psi[i].psi[mu];
	Psi[tid+dim]->psi[i].psi[mu].ptCV[ord] = Psi[tid+dim]->scfld->psi[i].psi[mu];
      }
    }
  }
#pragma omp barrier
#endif
}





#ifdef __RENORM_OMP__

void XiBuild1(ptGluon_fld& Umu, void** agg, int site, int col){

  ptSpinColor_fld **Psi;
  Psi = (ptSpinColor_fld**)agg;

  SpinColor Xi1[8],Xi2[8];
  ptSU3 *U;
  U = Umu.handle();
  int point_star;  

#pragma omp parallel num_threads(NTHR)
  {
    int *xx = new int[4];
    xx[0] =  1;
    xx[1] =  2;
    xx[2] =  3;
    xx[3] =  4;
    point_star = Umu.Z->get(xx);


    int tid = omp_get_thread_num();
    
    Psi[tid    ]->scfld->zeros();
    Psi[tid+dim]->scfld->zeros();
    
    memset(Psi[tid    ]->psi, 0, Psi[tid]->Z->Size*(sizeof(ptSpinColor)));
    memset(Psi[tid+dim]->psi, 0, Psi[tid]->Z->Size*(sizeof(ptSpinColor)));
    
    Psi[tid    ]->scfld->psi[site].psi[tid].whr[col].re = 1.0;
    Psi[tid+dim]->scfld->psi[site].psi[tid].whr[col].re = 1.0;
    
    Psi[tid    ]->M0invT();
    Psi[tid+dim]->M0inv();
    
    Psi[tid    ]->fftT(1);
    Psi[tid+dim]->fft(1);
    
    for( int i = 0 ; i < Umu.Z->Size; i++ ) {
      for (int mu = 0; mu < dim; mu++ ) {
	(Psi[tid    ]->psi[i]).psi[mu].ptCV[0] = ((Psi[tid    ]->scfld)->psi[i]).psi[mu];
	(Psi[tid+dim]->psi[i]).psi[mu].ptCV[0] = ((Psi[tid+dim]->scfld)->psi[i]).psi[mu];

      }
    }
 
#pragma omp barrier
    
    int point_c, link_c, slink_c, point_up, point_dn;
    SpinColor ZMom1, ZMom2;
 
  
    for( int jord = 1; jord <= params.ptord; jord++) {
    
      for( int c = 0 ; c < Umu.Z->Size; c++){
      
	point_c = Umu.Z->get(c);
      
	for( int kord = 0; kord < jord; kord++) {
	
	  for (int mu = 0; mu < dim; mu++){
	  
	    //trasposti
	    Psi[tid]->psi[point_c].psi[mu].ptCV[jord] += mcpt[jord-kord]*Psi[tid]->psi[point_c].psi[mu].ptCV[kord];
	  
	    //dritti
	    Psi[tid+dim]->psi[point_c].psi[mu].ptCV[jord] += mcpt[jord-kord]*Psi[tid+dim]->psi[point_c].psi[mu].ptCV[kord];
	  
	  }	
	
	  for( int mu = 0; mu < dim; mu++){
	  
	    link_c   = dim*(point_c)+mu;
	    slink_c  = dim*(Umu.Z->get(point_c, -1, mu))+mu;
	    point_up = Umu.Z->get(c, 1, mu);
	    point_dn = Umu.Z->get(c,-1, mu);
	  
	    for( int nu = 0; nu < dim; nu++ ){
	    
	      //trasposti
	      Xi1[tid].psi[nu] = Psi[tid]->psi[point_dn].psi[nu].ptCV[kord];
	      Xi2[tid].psi[nu] = Psi[tid]->psi[point_up].psi[nu].ptCV[kord];
	    
	      // dritti
	      Xi1[tid+dim].psi[nu] = Psi[tid+dim]->psi[point_dn].psi[nu].ptCV[kord];
	      Xi2[tid+dim].psi[nu] = Psi[tid+dim]->psi[point_up].psi[nu].ptCV[kord];
	    
	    }
	  
	    Xi1[tid    ] -= Xi1[tid    ].gmleft(mu)*gmT[mu];
	    Xi2[tid    ] += Xi2[tid    ].gmleft(mu)*gmT[mu];
	  
	    Xi1[tid+dim] += Xi1[tid+dim].gmleft(mu);
	    Xi2[tid+dim] -= Xi2[tid+dim].gmleft(mu);
	  
	    for( int nu = 0; nu < dim; nu++ ){
	    
	      // trasposti
	      Psi[tid    ]->psi[point_c].psi[nu].ptCV[jord] -= .5*tra(U[slink_c].ptU[jord-kord-1])*Xi1[tid    ].psi[nu];
	      Psi[tid    ]->psi[point_c].psi[nu].ptCV[jord] -= .5*(~(U[link_c].ptU[jord-kord-1]) )*Xi2[tid    ].psi[nu];
	    
	      // dritti
	      Psi[tid+dim]->psi[point_c].psi[nu].ptCV[jord] -= .5*dag(U[slink_c].ptU[jord-kord-1])*Xi1[tid+dim].psi[nu];
	      Psi[tid+dim]->psi[point_c].psi[nu].ptCV[jord] -= .5*( U[link_c].ptU[jord-kord-1] )*Xi2[tid+dim].psi[nu];
	    
	    } //nu
	  
	  } //mu

	}// kord

      }// sito
        
      for( int i = 0 ; i < Umu.Z->Size; i++ ) {
	for (int mu = 0; mu < dim; mu++ ) {
	  Psi[tid    ]->scfld->psi[i].psi[mu] = Psi[tid    ]->psi[i].psi[mu].ptCV[jord];
	  Psi[tid+dim]->scfld->psi[i].psi[mu] = Psi[tid+dim]->psi[i].psi[mu].ptCV[jord];
	}
      }
    
      // Subtract ZeroMom
      memset((void*)&ZMom1, 0, sizeof(SpinColor));
      memset((void*)&ZMom2, 0, sizeof(SpinColor));
      for( int i = 0 ; i < Umu.Z->Size; i++ ) {
	ZMom1 += Psi[tid+dim]->scfld->psi[i];
	ZMom2 += Psi[tid    ]->scfld->psi[i];
      }
      ZMom1 *= 1.0/(double)(Umu.Z->Size);
      ZMom2 *= 1.0/(double)(Umu.Z->Size);
      for( int i = 0 ; i < Umu.Z->Size; i++ ) {
	Psi[tid+dim]->scfld->psi[i] -= ZMom1;
	Psi[tid    ]->scfld->psi[i] -= ZMom2;
      }
    
      Psi[tid    ]->fftT(0);
      Psi[tid+dim]->fft(0);

      Psi[tid    ]->M0invT();
      Psi[tid+dim]->M0inv();
    
      Psi[tid    ]->fftT(1);
      Psi[tid+dim]->fft(1);
    
      for( int i = 0 ; i < Umu.Z->Size; i++ ) {
	for (int mu = 0; mu < dim; mu++ ) {
	  Psi[tid    ]->psi[i].psi[mu].ptCV[jord] = -1*Psi[tid    ]->scfld->psi[i].psi[mu];
	  Psi[tid+dim]->psi[i].psi[mu].ptCV[jord] = -1*Psi[tid+dim]->scfld->psi[i].psi[mu];
	
	}
      }
    } //jord
  }
#pragma omp barrier
  
  for( int jord = 0; jord <= params.ptord; jord++) {
    fft(Psi, 0, jord);
  }
  
}



#else



void XiBuild1(ptGluon_fld& Umu, void** agg, int site, int col){

  ptSpinColor_fld **Psi;
  Psi = (ptSpinColor_fld**)agg;

  SpinColor Xi1[8],Xi2[8];
  ptSU3 *U;
  U = Umu.handle();
  
  //  for( int tid = 0; tid < dim; tid++){
  int tid = 0;
  
  int *xx = new int[4];
  xx[0] = 1;
  xx[1] = 2;
  xx[2] = 3;
  xx[3] = 4;
  
  Psi[tid    ]->scfld->zeros();
  Psi[tid+dim]->scfld->zeros();
  
  memset(Psi[tid    ]->psi, 0, Psi[tid]->Z->Size*(sizeof(ptSpinColor)));
  memset(Psi[tid+dim]->psi, 0, Psi[tid]->Z->Size*(sizeof(ptSpinColor)));
  
  Psi[tid    ]->scfld->psi[site].psi[tid].whr[col].re = 1.0;
  Psi[tid+dim]->scfld->psi[site].psi[tid].whr[col].re = 1.0;
  
  Psi[tid    ]->M0invT();
  Psi[tid+dim]->M0inv();
  
  Psi[tid    ]->fftT(1);
  Psi[tid+dim]->fft(1);
  
  for( int i = 0 ; i < Umu.Z->Size; i++ ) {
    for (int mu = 0; mu < dim; mu++ ) {
      (Psi[tid    ]->psi[i]).psi[mu].ptCV[0] = ((Psi[tid    ]->scfld)->psi[i]).psi[mu];
      (Psi[tid+dim]->psi[i]).psi[mu].ptCV[0] = ((Psi[tid+dim]->scfld)->psi[i]).psi[mu];
      
    }
  }
  

  int point_c, link_c, slink_c, point_up, point_dn;
  SpinColor ZMom1, ZMom2;
  
  
  for( int jord = 1; jord <= params.ptord; jord++) {
    
    for( int c = 0 ; c < Umu.Z->Size; c++){
      
      point_c = Umu.Z->get(c);
      
      for( int kord = 0; kord < jord; kord++) {
	
	for (int mu = 0; mu < dim; mu++){
	  
	  //trasposti
	  Psi[tid]->psi[point_c].psi[mu].ptCV[jord] += mcpt[jord-kord]*Psi[tid]->psi[point_c].psi[mu].ptCV[kord];
	  
	  //dritti
	  Psi[tid+dim]->psi[point_c].psi[mu].ptCV[jord] += mcpt[jord-kord]*Psi[tid+dim]->psi[point_c].psi[mu].ptCV[kord];
	  
	}	
	
	for( int mu = 0; mu < dim; mu++){
	  
	  link_c   = dim*(point_c)+mu;
	  slink_c  = dim*(Umu.Z->get(point_c, -1, mu))+mu;
	  point_up = Umu.Z->get(c, 1, mu);
	  point_dn = Umu.Z->get(c,-1, mu);
	  
	  for( int nu = 0; nu < dim; nu++ ){
	    
	    //trasposti
	    Xi1[tid].psi[nu] = Psi[tid]->psi[point_dn].psi[nu].ptCV[kord];
	    Xi2[tid].psi[nu] = Psi[tid]->psi[point_up].psi[nu].ptCV[kord];
	    
	    // dritti
	    Xi1[tid+dim].psi[nu] = Psi[tid+dim]->psi[point_dn].psi[nu].ptCV[kord];
	    Xi2[tid+dim].psi[nu] = Psi[tid+dim]->psi[point_up].psi[nu].ptCV[kord];
	    
	  }
	  
	  Xi1[tid    ] -= Xi1[tid    ].gmleft(mu)*gmT[mu];
	  Xi2[tid    ] += Xi2[tid    ].gmleft(mu)*gmT[mu];
	  
	  Xi1[tid+dim] += Xi1[tid+dim].gmleft(mu);
	  Xi2[tid+dim] -= Xi2[tid+dim].gmleft(mu);
	  
	  for( int nu = 0; nu < dim; nu++ ){
	    
	    // trasposti
	    Psi[tid    ]->psi[point_c].psi[nu].ptCV[jord] -= .5*tra(U[slink_c].ptU[jord-kord-1])*Xi1[tid    ].psi[nu];
	    Psi[tid    ]->psi[point_c].psi[nu].ptCV[jord] -= .5*(~(U[link_c].ptU[jord-kord-1]) )*Xi2[tid    ].psi[nu];
	    
	    // dritti
	    Psi[tid+dim]->psi[point_c].psi[nu].ptCV[jord] -= .5*dag(U[slink_c].ptU[jord-kord-1])*Xi1[tid+dim].psi[nu];
	    Psi[tid+dim]->psi[point_c].psi[nu].ptCV[jord] -= .5*( U[link_c].ptU[jord-kord-1] )*Xi2[tid+dim].psi[nu];
	    
	  }
	  
	} //mu

	// if(tid == 0 && jord == 2 && point_c == site){
	//   (Psi[tid+dim]->psi[point_c]).psi[0].ptCV[jord].prout();
	//   (Psi[tid+dim]->psi[point_c]).psi[1].ptCV[jord].prout();
	//   (Psi[tid+dim]->psi[point_c]).psi[2].ptCV[jord].prout();
	//   (Psi[tid+dim]->psi[point_c]).psi[3].ptCV[jord].prout();
	//   cout << endl;
	// }
      
	
      }// kord

    }// sito
        
    for( int i = 0 ; i < Umu.Z->Size; i++ ) {
      for (int mu = 0; mu < dim; mu++ ) {
	Psi[tid    ]->scfld->psi[i].psi[mu] = Psi[tid    ]->psi[i].psi[mu].ptCV[jord];
	Psi[tid+dim]->scfld->psi[i].psi[mu] = Psi[tid+dim]->psi[i].psi[mu].ptCV[jord];
      }
    }
    
    // Subtract ZeroMom
    memset((void*)&ZMom1, 0, sizeof(SpinColor));
    memset((void*)&ZMom2, 0, sizeof(SpinColor));
    for( int i = 0 ; i < Umu.Z->Size; i++ ) {
      ZMom1 += Psi[tid+dim]->scfld->psi[i];
      ZMom2 += Psi[tid    ]->scfld->psi[i];
    }
    ZMom1 *= 1.0/(double)(Umu.Z->Size);
    ZMom2 *= 1.0/(double)(Umu.Z->Size);
    for( int i = 0 ; i < Umu.Z->Size; i++ ) {
      Psi[tid+dim]->scfld->psi[i] -= ZMom1;
      Psi[tid    ]->scfld->psi[i] -= ZMom2;
    }
    
    
    Psi[tid    ]->fftT(0);
    Psi[tid+dim]->fft(0);
    
    Psi[tid    ]->M0invT();
    Psi[tid+dim]->M0inv();
    
    Psi[tid    ]->fftT(1);
    Psi[tid+dim]->fft(1);
    
    for( int i = 0 ; i < Umu.Z->Size; i++ ) {
      for (int mu = 0; mu < dim; mu++ ) {
	Psi[tid    ]->psi[i].psi[mu].ptCV[jord] = -1*Psi[tid    ]->scfld->psi[i].psi[mu];
	Psi[tid+dim]->psi[i].psi[mu].ptCV[jord] = -1*Psi[tid+dim]->scfld->psi[i].psi[mu];
	
      }
    }
  
  } //jord

  
  for( int jord = 0; jord <= params.ptord; jord++) {
    fft(Psi, 0, jord);
  }
        
};



#endif


