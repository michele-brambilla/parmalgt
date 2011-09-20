/*
  Multithread, ma non ha molto senso farlo con piu' di 2 threads: il 32^4
  impiega ~45 min su singolo core

  Per l'utilizzo multithread definire la variabile _WILSON_LOOP_THREADS_
  Per misurare una configurazione di APE definire APE_CONFIG
*/

#define WHO_CARES_ABOUT_FERMIONS
#include <input.h>
#include <fstream>
#include <omp.h>

#ifdef NTHR
#undef NTHR
#endif
#define NTHR 2

nspt_params_t nspt_pars;
act_params_t  act_pars;

char *outn;
int *lato;
int half_len;
int eff_len;

double rho;

std::ofstream manyPlaq;

#ifndef _WILSON_LOOP_THREADS_
void measure(ptSU3*, ptSU3*, ptSU3*);
#else
void measure(ptSU3*, ptSU3*, ptSU3*, int);
#endif

void out(ptSU3*);
void out(ptSU3**);

void smearing(ptGluon_fld& Umu);


using namespace std;



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




int initialize(int argc, char **argv){
  
  FILE *fp;
  if( (fp = fopen("Loop.cfg","r") ) == NULL )
    {
      cout << "Impossibile leggere il file di configurazione.\nProgramma terminato." << endl;
      return 1;
    }
  fscanf(fp,"size %d %d %d %d\n", &act_pars.sz[0], &act_pars.sz[1], 
	 &act_pars.sz[2], &act_pars.sz[3]);
  get_val(fp, "conf", "%s" ,nspt_pars.confn );
  get_val(fp, "out" , "%s" ,outn            );
  fclose(fp);

  int opt;
  extern char *optarg;
  while ((opt = getopt(argc, argv, "l:c:o:r:h")) != -1) {
    switch (opt) {
    case 'l':
      act_pars.sz[0] = act_pars.sz[1] = act_pars.sz[2] = act_pars.sz[3] = atoi(optarg);
      break;
    case 'c':
      strcpy(nspt_pars.confn,optarg);
      break;
    case 'o':
      strcpy(outn,optarg);
      break;
    case 'h':
      cout << endl
	   << "-l\ttaglia (reticolo t==x==y==z);"  << endl
	   << "-c\tnome configurazione in input;"  << endl
	   << "-o\tnome file di output;"           << endl;
      exit(0);
      break;
    }
  }
  
  act_pars.iVol = act_pars.sz[0]*act_pars.sz[1]*act_pars.sz[2]*act_pars.sz[3];

  return 0;
}







int main(int argc, char** argv){

  lato = new int[dim];

  act_pars.sz     = new int[dim];
  nspt_pars.confn = new char[200];
  outn            = new char[200];

  Cplx *plaq = new Cplx[allocORD+1];

  if( initialize(argc, argv) ) exit(0);

#ifdef _WILSON_LOOP_THREADS_

#if (NTHR > 2)
#undef NTHR
#define NTHR 2
      cout << "Riduco il numero di threads a 2\n" << endl;
#endif

#endif

  cout << endl
       << "Taglia\t\t"       << act_pars.sz[0] << "\t" << act_pars.sz[1] << "\t" << act_pars.sz[2] << "\t" << act_pars.sz[3] << endl
       << "Size\t\t"         << act_pars.iVol     << endl
#ifdef _WILSON_LOOP_THREADS_
       << "Numero thread\t"  << NTHR     << endl
#endif
       << "Configurazione\t" << nspt_pars.confn << endl
       << "Output\t\t"       << outn            << endl;

  
  for(int mu = 0; mu < dim; ++mu)
    {
      lato[mu] = act_pars.sz[mu]/2;
    }
  
  
  latt LL(act_pars.sz);
  ptGluon_fld Umu(&LL);

#ifdef APE_CONFIG

  Umu.load_ape(nspt_pars.confn);

#else

  if( Umu.load(nspt_pars.confn,plaq) )
    {
      cout << "Impossibile leggere la configurazione." << endl;
      exit(-1);
    }

#endif


  manyPlaq.open(outn);
  manyPlaq.precision(15);

  cout << endl << "Placchetta:\t" << endl;
  for( int i1 = 2; i1 <= PTORD; i1+=2)
    {
      cout << plaq[i1].re << "\t";
    }
  cout << endl << endl;

  eff_len  = lato[0]+1;
  half_len = lato[0]/2;

  int curr, site_x, site_y;
#ifdef _WILSON_LOOP_THREADS_
  int tid;
#endif
  curr = 0;
  
  ptSU3 *appoggio0, *appoggio1, *tmp[dim*(dim-1)/2];
  appoggio0 = new ptSU3[eff_len*eff_len];
  appoggio1 = new ptSU3[eff_len*eff_len];

  for( int i = 0; i < dim*(dim-1)/2; i++)
    tmp[i] = new ptSU3[lato[0]*lato[1]];

  // Li inizializzo una volta per tutte
  appoggio0[0].flag = 1;
  appoggio1[0].flag = 1;
  for(int r = 1; r < eff_len; ++r)
    {
      appoggio1[r].flag = 1;
      appoggio0[eff_len*r].flag = 1;
    }

  int idx;

  for(curr = 0; curr < act_pars.iVol; curr++)
    {
      // mu-nu plane index
      idx = 0;
      for (int mu = 3; mu > 0; --mu)
	{
      
	  site_x = LL.L[curr][4];
      
	  // Riempio la colonna a nu=0 sulla componente 0
	  for(int r = 1; r < eff_len; ++r)
	    {
	      appoggio0[r] = appoggio0[r-1] * Umu.W[site_x].U[mu];
	      //	      appoggio1[r].flag = 1;
	      site_x = LL.L[site_x][5+mu];
	    }
	  
	  for(int nu = 0; nu < mu; ++nu)
	    {
	      // .. e fin qui..

	      // Riempio la colonna a mu=0
	      site_y = LL.L[curr][4];
	      for(int t = 1; t < eff_len; ++t)
		{
		  appoggio1[eff_len*t] = appoggio1[eff_len*(t-1)]*Umu.W[site_y].U[nu];
		  //		  appoggio0[eff_len*t].flag = 1;
		  site_y = LL.L[site_y][5+nu];
		}
	      // .. e un altro passo avanti..
		  
	      // riempio il bulk
		  
#ifndef _WILSON_LOOP_THREADS_
		  
	      // prima le colonne -> appoggio1
	      site_x = curr;
	      for(int r = 1; r < eff_len; ++r)
		{
		      
		  site_x = LL.L[site_x][5+mu];
		  site_y = site_x;
		  for(int t = 1; t < eff_len; ++t)
		    {
		      appoggio1[eff_len*t+r] = appoggio1[eff_len*(t-1)+r]*Umu.W[site_y].U[nu];
		      site_y = LL.L[site_y][5+nu];
		    } // end t
		      
		} // end r
		  
		  
		  // poi le righe -> appoggio0
	      site_y = curr;
	      for(int t = 1; t < eff_len; ++t)
		{
		  site_y = LL.L[site_y][5+nu];
		  site_x = site_y;
		  for(int r = 1; r < eff_len; ++r)
		    {
		      appoggio0[eff_len*t+r] = appoggio0[eff_len*t+r-1]*Umu.W[site_x].U[mu];
		      site_x = LL.L[site_x][5+mu];
		    } // end t
		      
		} // end r
		  
	      measure(appoggio0, appoggio1, tmp[idx]);
		  
#else
#pragma omp parallel private(site_x, site_y, tid) num_threads(NTHR)
	      {
		tid = omp_get_thread_num();
		    
		if( tid == 0)
		  {
		    // meta' dei thread si prende le colonne
			
		    site_x = curr;
		    for(int r = 1; r < eff_len; ++r)
		      {
			    
			site_x = LL.L[site_x][5+mu];
			site_y = site_x;
			for(int t = 1; t < eff_len; ++t)
			  {
			    appoggio1[eff_len*t+r] = appoggio1[eff_len*(t-1)+r]*Umu.W[site_y].U[nu];
			    site_y = LL.L[site_y][5+nu];
			  } // end t
			    
		      } // end r
		  }
		    
		// l'altra meta' le righe
		else
		  {
		    site_y = curr;
		    for(int t = 1; t < eff_len; ++t)
		      {
			site_y = LL.L[site_y][5+nu];
			site_x = site_y;
			for(int r = 1; r < eff_len; ++r)
			  {
			    appoggio0[eff_len*t+r] = appoggio0[eff_len*t+r-1]*Umu.W[site_x].U[mu];
			    site_x = LL.L[site_x][5+mu];
			  } // end t
			    
		      } // end r
		  }
		    
#pragma omp barrier
		    
		measure(appoggio0, appoggio1, tmp[idx], tid);
		    
	      }// end parallel
#endif	      
	      idx = 10*mu+nu;

	    } // end nu
	} // end mu
  
    } // end curr
  
  out(tmp);
  
  manyPlaq.close();
  return(0);
  
}

#ifndef _WILSON_LOOP_THREADS_


void measure(ptSU3* appoggio0, ptSU3* appoggio1, ptSU3 *tmp){

  for(int r = 1; r < eff_len; r++)
    {
      for(int t = 1; t < eff_len; t++)
	{
	  
	  tmp[(t-1)*lato[0]+r-1] += appoggio0[r]*appoggio1[t*eff_len+r]*dag(appoggio1[t*eff_len]*appoggio0[t*eff_len+r]);
	  
	} // end t
    } // end t
  
}


#else


void measure(ptSU3* appoggio0, ptSU3* appoggio1, ptSU3 *tmp, int tid){

  for(int r = 1+tid*half_len; r < (tid+1)*half_len+1; r++)
    {
      for(int t = 1; t < eff_len; t++)
	{
	  
	  tmp[(t-1)*lato[0]+r-1] += ( appoggio0[r]*appoggio1[t*eff_len+r]*
				      dag(appoggio1[t*eff_len]*appoggio0[t*eff_len+r]) );
	  
	} // end t
    } // end r
  
}

#endif




void out(ptSU3 *tmp){

  // normalizzazione : dim*(dim-1)/2 * NC *Vol
  double norm = 1.0/(double)(18.0*act_pars.iVol);

  for(int t = 1; t <= lato[0]; t++)
    {
      for(int r = 1; r <= lato[0]; r++)
	{
	  
	  manyPlaq << r << "\t"
		   << t << "\t";
	  for( int i1 = 1; i1 < PTORD; i1+=2)
	    {
	      manyPlaq << (tmp[(r-1)+(t-1)*lato[0]].ptU[i1].Tr()).re*norm << "\t";
	    }
	  manyPlaq << endl;
	    
	}
    }
  
}


void out(ptSU3 **tmp){

  // normalizzazione : dim*(dim-1)/2 * NC *Vol
  double norm = 1.0/(double)(18.0*act_pars.iVol);


  for(int id = 0; id < dim*(dim-1)/2; id++)
    {
      for(int t = 1; t <= lato[0]; t++)
	{
	  for(int r = 1; r <= lato[0]; r++)
	    {
	      
	      manyPlaq << id << "\t"
		       << r  << "\t"
		       << t  << "\t";
	      for( int i1 = 1; i1 < PTORD; i1+=2)
		{
		  manyPlaq << (tmp[id][(r-1)+(t-1)*lato[0]].ptU[i1].Tr()).re*norm << "\t";
		}
	      manyPlaq << endl;
	      
	    }
	}
    }
  
}


