#include<RenLib.h>

#include <sys/syscall.h>
#include <sys/types.h>
#include <unistd.h>
#include <omp.h>
#include <ctime>

#ifdef __MINIMAL_MPI__
#include "mpi.h"
#endif


using namespace std;

renorm_params_t params;


extern Cplx gmuval[15][4];
extern int gmuind[15][4];

static int get_val(FILE* fp, const char *str, const char* fmt,  void* val)
{
  char c[128];

  if(1!=fscanf(fp,"%s",c)){
    fprintf(stderr,"Error reading input file at %s\n",str);
    exit(1);
  }
    
  if(strcmp(str,c)!=0){
    fprintf(stderr,"Error reading input file expected %s found %s\n",str,c);
    exit(1);
  }

  if(1!=fscanf(fp,fmt,val)){
    fprintf(stderr,"Error reading input file at %s\n",str);
    fprintf(stderr,"Cannot read value format %s\n",fmt);
    exit(1);
  }
    
  return 0;
}





int get_params(int argc, char **argv){

  FILE *fin;
  if( (fin = fopen("Ren.cfg","r")) == NULL)
    {
      cout << "File di configurazione non trovato" << endl;
      return 1;
    }
  
  fscanf(fin, "taglia %d %d %d %d \n", &params.sz[0], &params.sz[1], &params.sz[2], &params.sz[3]);
  
  get_val(fin, "ordine", "%d",  &(params.ptord)  );
  get_val(fin, "conf"  , "%s",   (params.conf)   );
  get_val(fin, "out"   , "%s",   (params.out)    );
  get_val(fin, "prop"  , "%s",   (params.prfile) );
  get_val(fin, "alpha" , "%lf", &(params.alpha)  );
  
  fclose(fin);
  
  int opt;
  extern char *optarg;
  while ((opt = getopt(argc, argv, "l:o:d:c:a:p:f:r:h")) != -1) {
    switch (opt) {
    case 'l':
      params.sz[0] = params.sz[1] = params.sz[2] =
	params.sz[3] = atoi(optarg);
      break;
    case 'o':
      params.ptord = atoi(optarg);
      break;
    case 'c':
      strcpy(params.conf,optarg);
      break;
    case 'f':
      strcpy(params.out,optarg);
      break;
    case 'p':
      strcpy(params.prfile,optarg);
      break;
    case 'a':
      params.alpha = atof(optarg);
      break;
    case 'h':
      cout << endl
	   << "-l\t\ttaglia (reticolo t==x==y==z)\n"
	   << "-o\t\tordine perturbativo\n"
	   << "-c\t\tnome configurazione\n"
	   << "-f\t\tnome file output\n"
	   << "-p\t\tnome output propagatore\n"
	   << "-a\t\talpha\n"
	   << "-t\t\tnthr\n"
	   << endl;
      exit(0);
      break;
    }
  }
  
  params.iVol = params.sz[0] * params.sz[1] * params.sz[2] * params.sz[3];

  
  uint i = 0;
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
  if( i == strlen(params.conf) )
    {
      std::cout << std::endl << "Wrong gauge action." << std::endl << std::endl;
      return 1;
    }



#ifdef __MINIMAL_MPI__
  int rank;
  char *tmp_conf  = new char[110];
  char *tmp_conf1 = new char[110];
  int len_str = strlen(params.conf);
  {
    int rc = MPI_Init(&argc,&argv);
    if (rc != MPI_SUCCESS) {
      printf ("Error starting MPI program. Terminating.\n");
      MPI_Abort(MPI_COMM_WORLD, rc);
    }
    
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
  }

  while(len_str > 0 )
    {
      if(strncmp(params.conf+len_str,".",1) == 0) break;   
      len_str--;
    }
  strncpy(tmp_conf,params.conf,len_str);
  //  printf("%s.r%d%s\n",tmp_conf,rank,params.conf+len_str);
  sprintf(tmp_conf1,"%s.r%d%s",tmp_conf,rank,params.conf+len_str);
  strcpy(params.conf,tmp_conf1);
  
  delete [] tmp_conf;
  delete [] tmp_conf1;
#endif


  strcat(params.prfile,"_");
  strcat(params.prfile,params.conf+i);

  strcat(params.out,"_");
  strcat(params.out,params.conf+i);

#ifdef __TRANSP_PROPAG__
  strncpy(params.prfileT,params.prfile,strlen(params.prfile)-4);
  strcat(params.prfileT,".tra");
#endif

  strcat(params.mom,"running.log.");
  strcat(params.mom,params.conf+i);
  
  

  // Verifica taglia reticolo
  char* lsize = new char[20];
  sprintf(lsize,"%dx%dx%dx%d",
	  params.sz[0],params.sz[1],params.sz[2],params.sz[3]);
  cout << lsize << endl;

  i = 0;

  if(params.sz[0] < 10)
    {
      while(i < strlen(params.conf))
	{
	  if( !strncmp(lsize,params.conf+i,7) ) break;
	  i++;
	}
    }
  else
    if(params.sz[0] < 100)
      {
	while(i < strlen(params.conf))
	  {
	    if( !strncmp(lsize,params.conf+i,11) ) break;
	    i++;
	  }
      }
  
  if( i == strlen(params.conf) )
    {
      cout << endl << "Wrong lattice size." << endl << endl;
      return 1;
    }
  delete [] lsize;

  
  cout << "Parametri:\n" 
       << "taglia " <<params.sz[0]<<" "<< params.sz[1] <<" "<< params.sz[2] <<" "<< params.sz[3] <<"\n"
       << "ordine "         << params.ptord   << "\n"
       << "configurazione " << params.conf    << "\n"
       << "output "         << params.out     << "\n"
       << "propagatore "    << params.prfile  << "\n"
#ifdef __TRANSP_PROPAG__
       << "propagatoreTra " << params.prfileT << "\n"
#endif
       << "alpha "          << params.alpha   << "\n";  

  params.momout << "Configurazione: " << params.conf << endl;
  params.momout << "File Output Corrente: " << params.out << endl;
  params.momout << "File Propagatore: " << params.prfile << endl;
#ifdef __TRANSP_PROPAG__
  params.momout << "File PropagatoreTrasposto: " << params.prfileT << endl << endl;
#endif
  
  time_t tempo;
  time ( &tempo );
  params.momout << "Inizio run: " << ctime(&tempo) << endl << endl;
  params.momout << "Momenti eseguiti" << endl;

  return 0;
}




// Misura a configurazione congelata la placchetta
int plaquette_measure(ptGluon_fld& Umu, Cplx* w2){
#ifndef __RENORM_OMP__
  
  Cplx* app = new Cplx[PTORD+1];
  Cplx* w = new Cplx[PTORD+1];
  ptSU3 W1, W;
  int link_c;

  for(int i = 0; i < params.iVol;i++){
    //    W1.zero();				
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
    w[i1] /= (params.iVol*72);
  }					

  delete [] app;

#else

  Cplx* app  = new Cplx[PTORD+1];
  Cplx* w = new Cplx[PTORD+1];
  Cplx** ww   = new Cplx*[NTHR];
  Cplx** ww1  = new Cplx*[NTHR];
  ptSU3 Ww1[NTHR], Ww[NTHR];
  int tid;

  for(int nt = 0; nt < NTHR; nt++){
    ww[nt]  = new Cplx[PTORD+1];
    ww1[nt] = new Cplx[PTORD+1];
  }
  
#pragma omp parallel private(tid) num_threads(NTHR)
  {
    tid = omp_get_thread_num();

#pragma omp for
    for(int site_x = 0; site_x < params.iVol;site_x++){
      Ww1[tid].zero();
      for(int mu = 0; mu < dim; mu++){
 	Ww[tid].zero();
 	for(int nu = 0; nu < dim; nu++){
 	  if(nu != mu ){
 	    Ww[tid] += Umu.staple(site_x, mu, nu);
 	  }
 	}
 	Ww1[tid] = Umu.W[site_x].U[mu]*Ww[tid];
 	Ww1[tid].Tr(ww1[tid]);
	
 	for(int i1 = 0; i1 <= PTORD; i1++){
 	  ww[tid][i1] += ww1[tid][i1];
 	}
	
      }// altro link, nello stesso sito	
    } // fine siti
    
  } // end parallel
#pragma omp barrier

  for(int i1 = 0; i1 <= PTORD; i1++){
    w[i1] = 0;
    for(int nt = 0; nt < NTHR; nt++){
      //      std::cout << ww[nt][i1].re << "\t";
      w[i1] += ww[nt][i1];
    }
    //    std::cout << std::endl;
    w[i1] /= (params.iVol*72);		
  }					
  std::cout << std::endl;

  delete [] app;

  delete [] ww;
  delete [] ww1;

#endif

  std::cout << params.iVol << std::endl;
  for(int i1 = 0; i1 <= PTORD;i1++){
    w[i1].prout();
    std::cout << "\t\t";
    w2[i1].prout();
    std::cout << std::endl;
  }

  for(int i1 = 0; i1 <= PTORD;i1++){
    if( (w[i1]-w2[i1]).mod() > 1e-7 ) return 1;
  }
  return 0; 

}








int main(int argc, char **argv){



  params.sz = new int(dim);
  params.conf    = new char[100];
  params.out     = new char[100];
  params.prfile  = new char[100];
#ifdef __TRANSP_PROPAG__
  params.prfileT = new char[100];
#endif
  params.mom     = new char[100];

  if(get_params(argc, argv)){
    //    exit(-1);
  }


  latt LL(params.sz);
  LL.p_init();

  ptGluon_fld Umu(&LL);

  cout.precision(10);

  Cplx* w2 = new Cplx[PTORD+1];
  if( carica_campo(Umu,w2) )
    {
      exit(-1);
    }

  if (plaquette_measure(Umu, w2) )
    {
      cout << "Errore nella lettura della configurazione." << endl;
      return 1;
    }
  delete [] w2;

  
  int stat, cnt;
  
  params.momout.open(params.mom);
    
  if(gauge_fixing(Umu, cnt)){
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


  Cplx propag[dim][dim][allocORD+1], curr[16][dim][dim][NC][allocORD+1];
#ifdef __TRANSP_PROPAG__
  Cplx propagT[dim][dim][allocORD+1];
#endif

  ptSpinColor gmp[16][4];


  int c;

  FILE *fout, *pout;
  fout = fopen(params.out,"wb");
  pout = fopen(params.prfile,"wb");
#ifdef __TRANSP_PROPAG__
  FILE *tout;
  tout = fopen(params.prfileT,"wb");
#endif
  


  FILE *fmom;
  if( (fmom = fopen("momREF.txt","r") ) == NULL){
    cout << "Errore, Momenti di Riferimento Non Trovati!" << endl;
    exit(1);
  }  
  int *xx       = new int[4];

  // ofstream of0, of1, of2, of3, of4, of5, of6, of7;
  // of0.open("Pmu0s.txt");
  // of1.open("Pmu1s.txt");
  // of2.open("Pmu2s.txt");
  // of3.open("Pmu3s.txt");
  // of4.open("Pmu4s.txt");
  // of5.open("Pmu5s.txt");
  // of6.open("Pmu6s.txt");
  // of7.open("Pmu7s.txt");





  while( ( stat=fscanf(fmom,"%d %d %d %d\n", &xx[0], &xx[1], &xx[2], &xx[3])) >0){ ;
    
    c = LL.get(xx);
    cout << "x = \t"    <<  xx[0] << "\t"
  	 << "y = \t"    <<  xx[1] << "\t"
  	 << "z = \t"    <<  xx[2] << "\t"
  	 << "t = \t"    <<  xx[3] << "\t"
  	 << "sito = \t" <<  c     << "\t"
  	 << "stato = "  <<  stat  << endl;
    
    
    bzero(propag, dim*dim*(allocORD+1)*sizeof(Cplx));
    bzero(curr  , 16*dim*dim*NC*(allocORD+1)*sizeof(Cplx));
    // memset(propag,0 ,sizeof(propag));
    // memset(curr,0 ,sizeof(curr));
#ifdef __TRANSP_PROPAG__
    //    memset(propagT,0 ,sizeof(propagT));
    bzero(propagT, dim*dim*(allocORD+1)*sizeof(Cplx));
#endif
    
    
    for(int a = 0; a < NC; a++){
      
      XiBuild1(Umu, (void**)Pmu, c, a);


      // Pmu[0]->fout(of0,2);
      // Pmu[1]->fout(of1,2);
      // Pmu[2]->fout(of2,2);
      // Pmu[3]->fout(of3,2);
      // Pmu[4]->fout(of4,2);
      // Pmu[5]->fout(of5,2);
      // Pmu[6]->fout(of6,2);
      // Pmu[7]->fout(of7,2);

#ifdef __RENORM_OMP__
#pragma omp barrier
#pragma omp parallel for num_threads(NTHR)  
#endif      
      for(int mu = 0; mu < dim; mu++) {
	for(int o = 0 ; o <= params.ptord; o++){
#ifdef __TRANSP_PROPAG__	  
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
    
    
    
#ifdef __RENORM_OMP__
#pragma omp barrier
#pragma omp parallel for num_threads(NTHR)
#endif      
      for(int cid = 0; cid < 16; cid++)
	{
	
	  for( int d = 0; d < LL.Size; d++ ) {
	  
	    // Apply gammas
	    for(int i = 0; i < dim; i++)  
	      {
	    	for(int nu = 0; nu < dim; nu++)
	    	  {
		    
	    	    if( cid > 0 )
	    	      {
	    		gmp[cid][i].psi[nu] = gmuval[cid-1][nu]*(Pmu[i+dim]->psi[d].psi[gmuind[cid-1][nu]]);
	    	      }
	    	    else
	    	      {
	    		gmp[cid][i].psi[nu] = (Pmu[i+dim]->psi[d].psi[nu]);
	    	      }
		    
	    	  }
	      }
	    
	    // In gmleft notation
	    // if(cid == 0)
	    //   {
	    // 	for(int i = 0; i < dim; i++)  
	    // 	  {
	    // 	    gmp[cid][i] = Pmu[i+dim]->psi[d];
	    // 	  }
	    //   }
	    // if(cid > 0 && cid < 5)
	    //   {
	    // 	for(int i = 0; i < dim; i++)  
	    // 	  {
	    // 	    gmp[cid][i] = (Pmu[i+dim]->psi[d]).gmleft(cid-1);
	    // 	  }
	    //   }
	    // if(cid == 5)
	    //   {
	    // 	for(int i = 0; i < dim; i++)  
	    // 	  {
	    // 	    gmp[cid][i] = (Pmu[i+dim]->psi[d]).gmleft(5);
	    // 	  }
	    //   }
	    // if(cid > 5 && cid < 10)
	    //   {
	    // 	for(int i = 0; i < dim; i++)  
	    // 	  {
	    // 	    gmp[cid][i] = ((Pmu[i+dim]->psi[d]).gmleft(5)).gmleft(cid-6);
	    // 	  }
	    //   }

	    // if(cid == 10)
	    //   {
	    // 	for(int i = 0; i < dim; i++)  
	    // 	  {
	    // 	    gmp[cid][i] = ((Pmu[i+dim]->psi[d]).gmleft(1)).gmleft(0);
	    // 	  }
	    //   }
	    // if(cid == 11)
	    //   {
	    // 	for(int i = 0; i < dim; i++)  
	    // 	  {
	    // 	    gmp[cid][i] = ((Pmu[i+dim]->psi[d]).gmleft(2)).gmleft(0);
	    // 	  }
	    //   }
	    // if(cid == 12)
	    //   {
	    // 	for(int i = 0; i < dim; i++)  
	    // 	  {
	    // 	    gmp[cid][i] = ((Pmu[i+dim]->psi[d]).gmleft(3)).gmleft(0);
	    // 	  }
	    //   }
	    // if(cid == 13)
	    //   {
	    // 	for(int i = 0; i < dim; i++)  
	    // 	  {
	    // 	    gmp[cid][i] = ((Pmu[i+dim]->psi[d]).gmleft(2)).gmleft(1);
	    // 	  }
	    //   }
	    // if(cid == 14)
	    //   {
	    // 	for(int i = 0; i < dim; i++)  
	    // 	  {
	    // 	    gmp[cid][i] = ((Pmu[i+dim]->psi[d]).gmleft(3)).gmleft(1);
	    // 	  }
	    //   }
	    // if(cid == 15)
	    //   {
	    // 	for(int i = 0; i < dim; i++)  
	    // 	  {
	    // 	    gmp[cid][i] = ((Pmu[i+dim]->psi[d]).gmleft(3)).gmleft(2);
	    // 	  }
	    //   }


	  
	    for(int mu = 0; mu < dim; mu++) {
	      for(int nu = 0; nu < dim; nu++) {
	      
		for( int iO = 0; iO <=params.ptord; iO++ ){
		  for( int jO = 0; jO <= iO; jO++ ){
		  
		    for( int rho = 0; rho < dim; rho++) {
		      for ( int icol = 0; icol < NC; icol++ ) {
		      
			curr[cid][mu][nu][a][iO] += 
			  Pmu[mu]->psi[d].psi[rho].ptCV[jO].whr[icol] * 
			  gmp[cid][nu].psi[rho].ptCV[iO-jO].whr[icol] ;
			
		      } // color index saturation
		    } // dirac index saturation
		  
		  } // jord
		} // iord
		
	      } // nu
	    } // mu
	  
	  } // momentum index saturation
	
	} // current id (cid)
      
    } // color
  
  
  
    fwrite(curr,    sizeof(curr),    1, fout);
    fwrite(propag,  sizeof(propag),  1, pout);
#ifdef __TRANSP_PROPAG__
    fwrite(propagT, sizeof(propagT), 1, tout);
#endif
    params.momout << xx[0] << "\t"
		  << xx[1] << "\t"
		  << xx[2] << "\t"
		  << xx[3] << "\n";
    params.momout.flush();
  
  }
  
  fclose(fmom);
  
  time_t tempo;
  time ( &tempo );
  params.momout << "Fine run: " << ctime(&tempo) << endl << endl;
  
  
  fclose(fout);
  fclose(pout);
#ifdef __TRANSP_PROPAG__
  fclose(tout);
#endif
  params.momout.close();
  return 0;
}

