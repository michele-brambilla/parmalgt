#include "RenLib.h"

#include <sys/syscall.h>
#include <sys/types.h>
#include <unistd.h>
#include <omp.h>

using namespace std;

renorm_params_t params;

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
  if( (fin = fopen("GF.cfg","r")) == NULL)
    {
      cout << "File di configurazione non trovato" << endl;
      return 1;
    }
  
  fscanf(fin, "taglia %d %d %d %d \n", &params.sz[0], &params.sz[1], &params.sz[2], &params.sz[3]);
  
  get_val(fin, "ordine", "%d",  &(params.ptord)  );
  get_val(fin, "conf"  , "%s",   (params.conf)   );
  get_val(fin, "out"   , "%s",   (params.out)    );
  get_val(fin, "alpha" , "%lf", &(params.alpha)  );
  
  fclose(fin);
  
  int opt;
  extern char *optarg;
  while ((opt = getopt(argc, argv, "l:o:d:c:a:f:h")) != -1) {
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
    case 'a':
      params.alpha = atof(optarg);
      break;
    case 'h':
      cout << endl
	   << "-l\t\ttaglia (reticolo t==x==y==z)\n"
	   << "-o\t\tordine perturbativo\n"
	   << "-c\t\tnome configurazione\n"
	   << "-f\t\tnome file output\n"
	   << "-a\t\talpha\n"
	   << endl;
      exit(0);
      break;
    }
  }
  
  uint i = 0;

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
       << "alpha "          << params.alpha   << "\n"
       << endl;

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

  for(int i1 = 0; i1 <= PTORD;i1++){
    if( (w[i1]-w2[i1]).mod() > 1e-7 ) return 1;
  }
  return 0; 

}



int plaquette_measure(ptGluon_fld& Umu, Cplx* w2, Cplx* w){

  bzero(w,allocORD*sizeof(Cplx));
#ifndef __RENORM_OMP__
  
  Cplx* app = new Cplx[PTORD+1];
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

  for(int i1 = 0; i1 <= PTORD;i1++){
    if( (w[i1]-w2[i1]).mod() > 1e-7 ) return 1;
  }
  return 0; 

}








int main(int argc, char **argv){

  params.sz = new int(dim);
  params.conf    = new char[100];
  params.out     = new char[100];

  if(get_params(argc, argv)){
    exit(-1);
  }

  latt LL(params.sz);
  LL.p_init();

  ptGluon_fld Umu(&LL);

  cout.precision(10);

  Cplx* w2 = new Cplx[allocORD+1];
  Cplx* w  = new Cplx[allocORD+1];

  if( carica_campo(Umu,w2) )
    {
      exit(-1);
    }

  if (plaquette_measure(Umu, w2) )
    {
      cout << "Errore nella lettura della configurazione." << endl;
      return 1;
    }

  int cnt = 0;
  if(gauge_fixing(Umu, cnt)){
    cout << "Attenzione, il gauge fixing non converge" << endl;
    exit(-1);
  }


  if (plaquette_measure(Umu, w2, w) )
    {
      cout << "Errore nella gauge fixing." << endl;
      exit(-1);
    }
  delete [] w2;


  if( Umu.save(params.conf, w) ){
     cout << "Errore in salvataggio." << endl;
     exit(-1);
   }
  delete [] w;

  return 0;
}

