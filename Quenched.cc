#include<fstream>
#include<string>
#include <fenv.h>

#include"nspt.h"
#include "MassC.h"

using namespace std;

#ifdef __MINIMAL_MPI__
#include "mpi.h"
#endif

Cplx* w       = new Cplx[PTORD+1];
Cplx* w1      = new Cplx[PTORD+1];
Cplx* w2      = new Cplx[PTORD+1];
double* norm  = new double[PTORD];
double* norm1 = new double[PTORD];
Cplx** trU;
ptSU3* U;

nspt_params_t   nspt_pars;
act_params_t    act_pars;
thread_params_t thr_pars[NTHR];

int nInt = 0;

MyRand Rand[NTHR];


std::ofstream plaqfile, normfile, logfile, trufile;

#ifdef __WILLOOP_2x2__
std::ofstream loop2x2;
#endif

#ifdef __TIMING__
PRlgtTime Time;
#endif


#ifdef __K_MOM_ANALYSIS__
extern ptSU3_fld* Wgauge;
extern fftw_plan *planFA;

extern fftw_plan *planGluon;

extern double *knorm[allocORD+1];
FILE *FPkn;
#endif

// Lettura sicura dei parametri, verifica che la stringa che
// precede il valore del parametro corrisponda a quanto atteso

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




// Riceve i parametri da file e/o da linea di comando,
// apre i file necessari

int initialize(int argc, char** argv){

  long unsigned seed[NTHR];
  act_pars.sz       = new int[dim];
  nspt_pars.plaqn   = new char[100];
  nspt_pars.confn   = new char[100];
  nspt_pars.damon   = new char[100];
  nspt_pars.normn   = new char[100];
  nspt_pars.logn    = new char[100];
  nspt_pars.trun    = new char[100];
#ifdef __WILLOOP_2x2__
  nspt_pars.name2x2 = new char[100];
#endif

  trU = new Cplx* [dim];
  for( int i = 0; i < dim; i++)
    {
      trU[i] = new Cplx[PTORD];
    }


  // Legge i parametri da file. Se non ci riesce, 
  // ritorna un segnale d'errore
  FILE *fp;
  if( (fp = fopen("Quench.cfg","r")) == NULL){
    std::cout << "Errore: impossibile leggere il file di configurazione." 
	      << std::endl;
    return 1;
  }

  fscanf( fp,"taglia %d %d %d %d\n", &act_pars.sz[0], &act_pars.sz[1], 
	  &act_pars.sz[2], &act_pars.sz[3] );

  get_val(fp, "SWEEP",       "%d" ,&(nspt_pars.Sweep) );
  get_val(fp, "BEAT",        "%d" ,&(nspt_pars.Beat)  );
  get_val(fp, "tau_g",       "%lf",&(act_pars.tau_g)  );
  get_val(fp, "alpha",       "%lf",&(act_pars.alpha)  );
  get_val(fp, "init_status", "%d" ,&(nspt_pars.Init)  );
  get_val(fp, "plaq_out",    "%s" ,nspt_pars.plaqn    );
  get_val(fp, "last_conf",   "%s" ,nspt_pars.confn    );
  get_val(fp, "PTORD",       "%d" ,&(PTORD)           );
  for( int tid = 0; tid < NTHR; tid++)
    {
      get_val(fp, "seed",        "%ld",&seed[tid]     );
    }

  fclose(fp);

  // File damocle
  strcpy(nspt_pars.damon,"damocle.dag");


  // Sovrascrive le impostazioni da file con eventuali
  // parametri da riga di comando

  int opt;
  extern char *optarg;
  while ((opt = getopt(argc, argv, "l:s:b:d:r:g:a:t:p:c:n:h")) != -1) {
    switch (opt) {
    case 'l':
      act_pars.sz[0] = act_pars.sz[1] = act_pars.sz[2] = act_pars.sz[3] = atoi(optarg);
      break;
    case 's':
      nspt_pars.Sweep = atoi(optarg);
      break;
    case 'b':
      nspt_pars.Beat = atoi(optarg);
      break;
    case 'r':
#ifndef __PARALLEL_OMP__
      seed[0] = atol(optarg);
#else
      srand(atol(optarg));
      for(int thr = 0; thr < NTHR; thr++){
	seed[thr] = rand();
      }      
#endif
      break;
    case 'n':
      PTORD = atoi(optarg);
      break;
    case 'g':
      act_pars.tau_g = atof(optarg);
      break;
    case 'a':
      act_pars.alpha = -atof(optarg);
      break;
    case 't':
      nspt_pars.Init = atoi(optarg);
      break;
    case 'p':
      strcpy(nspt_pars.plaqn,optarg);
      break;
    case 'c':
      strcpy(nspt_pars.confn,optarg);
      break;
    case 'd':
      strcpy(nspt_pars.damon,optarg);
      break;
    case 'h':
      printf( "\n" );
      printf( "-l\t\ttaglia (reticolo t==x==y==z)\n" );
      printf( "-s\t\tnumero iterazioni (SWEEP)\n" );
      printf( "-b\t\titerazioni tra letture damocle (BEAT)\n" );
      printf( "-r\t\tseed generatore random\n" );
      printf( "-g\t\ttau gluoni\n" );
      printf( "-a\t\talpha\n" );
      printf( "-t\t\tstato iniziale: 1)da configurazione; 2)freddo\n" );
      printf( "-p\t\tnome file placchetta\n" );
      printf( "-c\t\tnome configurazione salvata\n" );
      printf( "-d\t\tnome file damocle\n" );
      printf( "-n\t\tordine pt\n" );
      printf( "\n" );
      exit(0);
      break;
    }
  }


  {
    FILE* df;
    if( ( df = fopen(nspt_pars.damon,"r") ) == NULL)
      {
	std::cout << "Errore, impossibile trovare file damocle." << std::endl;
	fclose(df);
	exit(-1);
      }
    fclose(df);
  }



  // Utilizza MPI per una farm
#ifdef __MINIMAL_MPI__
  int rank;
  int rc = MPI_Init(&argc,&argv);
  if (rc != MPI_SUCCESS) {
    printf ("Error starting MPI program. Terminating.\n");
    MPI_Abort(MPI_COMM_WORLD, rc);
  }

  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
#endif

  // Inizializza il generatore random
#ifdef __MINIMAL_MPI__

  time(NULL);
  for(int thr = 0; thr < NTHR; thr++){
    Rand[thr].init( seed[thr]+clock() );
  }

#else

  for(int thr = 0; thr < NTHR; thr++){
    Rand[thr].init(seed[thr]);
  }

#endif


  // Scelta dei coefficienti per le azioni improved
#if GAUGE_ACTION == WIL
  act_pars.c1 = 0.0;
#elif GAUGE_ACTION == TLSYM
  act_pars.c1 = -1.0/12.0;  
#elif GAUGE_ACTION == IWA
  act_pars.c1 = -0.331;  
#elif GAUGE_ACTION == DBW2
  act_pars.c1 = -1.4088;  
#endif
  act_pars.c0 = 1. - 8. * act_pars.c1;



  // Compone il nome dei files
  char *type_g = new char[3];
#if GAUGE_ACTION == WIL
  strcpy( type_g, "wil" );
#elif GAUGE_ACTION == TLSYM
  strcpy( type_g, "tls" );
#elif GAUGE_ACTION == IWA
  strcpy( type_g, "iwa" );
#elif GAUGE_ACTION == DBW2
  strcpy( type_g, "dbw" );
#endif

#ifdef __MINIMAL_MPI__
  
#ifdef __WILLOOP_2x2__
  sprintf(nspt_pars.name2x2,"%s2x2_%s_nf0_%dx%dx%dx%d_a%2.2f_tg%3.3f_o%d.r%d.txt",
	  nspt_pars.plaqn, type_g, act_pars.sz[0], act_pars.sz[1], 
	  act_pars.sz[2], act_pars.sz[3], fabs(act_pars.alpha), act_pars.tau_g,
	  PTORD,rank);
#endif

  sprintf(nspt_pars.plaqn,"%s_%s_nf0_%dx%dx%dx%d_a%2.2f_tg%3.3f_o%d.r%d.txt",
	  nspt_pars.plaqn, type_g, act_pars.sz[0], act_pars.sz[1], 
	  act_pars.sz[2], act_pars.sz[3], fabs(act_pars.alpha), act_pars.tau_g,
	  PTORD,rank);
#ifndef APE_CONFIG
  sprintf(nspt_pars.confn,"%s%s_nf0_%dx%dx%dx%d_a%2.2f_tg%3.3f_o%d.r%d.dat",
	  nspt_pars.confn, type_g, act_pars.sz[0], act_pars.sz[1], 
	  act_pars.sz[2],act_pars.sz[3],fabs(act_pars.alpha), act_pars.tau_g,
	  allocORD,rank);
#endif
  sprintf(nspt_pars.logn,"%s_nf0_%dx%dx%dx%d_a%2.2f_tg%3.3f_o%d.r%d.log",
	  type_g, act_pars.sz[0], act_pars.sz[1], 
	  act_pars.sz[2],act_pars.sz[3],fabs(act_pars.alpha), act_pars.tau_g,
	  PTORD,rank);
  sprintf(nspt_pars.normn,"%s_nf0_%dx%dx%dx%d_a%2.2f_tg%3.3f_o%d.r%d.nor",
	  type_g, act_pars.sz[0], act_pars.sz[1], 
	  act_pars.sz[2],act_pars.sz[3],fabs(act_pars.alpha), act_pars.tau_g, 
	  PTORD,rank);
  sprintf(nspt_pars.trun,"%s_nf0_%dx%dx%dx%d_a%2.2f_tg%3.3f_o%d.r%d.trU",
	  type_g, act_pars.sz[0], act_pars.sz[1], 
	  act_pars.sz[2],act_pars.sz[3],fabs(act_pars.alpha), act_pars.tau_g,
	  PTORD,rank);

#else

#ifdef __WILLOOP_2x2__
  sprintf(nspt_pars.name2x2,"%s2x2_%s_nf0_%dx%dx%dx%d_a%2.2f_tg%3.3f_o%d.txt",
	  nspt_pars.plaqn, type_g, act_pars.sz[0], act_pars.sz[1], 
	  act_pars.sz[2], act_pars.sz[3], fabs(act_pars.alpha), act_pars.tau_g,
	  PTORD);
#endif

  sprintf(nspt_pars.plaqn,"%s_%s_nf0_%dx%dx%dx%d_a%2.2f_tg%3.3f_o%d.txt",
	  nspt_pars.plaqn, type_g, act_pars.sz[0], act_pars.sz[1], 
	  act_pars.sz[2], act_pars.sz[3], fabs(act_pars.alpha), act_pars.tau_g,
	  PTORD);
#ifndef APE_CONFIG
  sprintf(nspt_pars.confn,"%s%s_nf0_%dx%dx%dx%d_a%2.2f_tg%3.3f_o%d.dat",
	  nspt_pars.confn, type_g, act_pars.sz[0], act_pars.sz[1], 
	  act_pars.sz[2],act_pars.sz[3],fabs(act_pars.alpha), act_pars.tau_g,
	  allocORD);
#endif
  sprintf(nspt_pars.logn,"%s_nf0_%dx%dx%dx%d_a%2.2f_tg%3.3f_o%d.log",
	  type_g, act_pars.sz[0], act_pars.sz[1], 
	  act_pars.sz[2],act_pars.sz[3],fabs(act_pars.alpha), act_pars.tau_g,
	  PTORD);
  sprintf(nspt_pars.normn,"%s_nf0_%dx%dx%dx%d_a%2.2f_tg%3.3f_o%d.nor",
	  type_g, act_pars.sz[0], act_pars.sz[1], 
	  act_pars.sz[2],act_pars.sz[3],fabs(act_pars.alpha), act_pars.tau_g,
	  PTORD);
  sprintf(nspt_pars.trun,"%s_nf0_%dx%dx%dx%d_a%2.2f_tg%3.3f_o%d.trU",
	  type_g, act_pars.sz[0], act_pars.sz[1], 
	  act_pars.sz[2],act_pars.sz[3],fabs(act_pars.alpha), act_pars.tau_g,
	  PTORD);

#endif

  // Apre i file su cui scrivere i dati
  if(!plaqfile.is_open() ) plaqfile.open(nspt_pars.plaqn, std::ios_base::app);
  if(!normfile.is_open() ) normfile.open(nspt_pars.normn);
  if(!logfile.is_open()  ) logfile.open(nspt_pars.logn);
  if(!trufile.is_open()  ) trufile.open(nspt_pars.trun);
#ifdef __WILLOOP_2x2__
  if(!loop2x2.is_open() ) loop2x2.open(nspt_pars.name2x2, std::ios_base::app);
  loop2x2.precision(20);
#endif

  // Dati di inizializzazione su su logfile
  logfile << "Taglia\t\t"       << act_pars.sz[0]  << " " << act_pars.sz[1] << " " << act_pars.sz[2] << " " << act_pars.sz[3] << std::endl
	  << "Sweep\t\t"        << nspt_pars.Sweep   << std::endl
	  << "Beat\t\t"         << nspt_pars.Beat    << std::endl
	  << "Tau_g\t\t"        << act_pars.tau_g    << std::endl
	  << "Alpha\t\t"        << act_pars.alpha    << std::endl
	  << "Init status\t"    << nspt_pars.Init    << std::endl
	  << "Ordine\t\t"       << PTORD             << std::endl
	  <<                                            std::endl
#ifdef __PARALLEL_OMP__
	  << "Numero Threads\t" << NTHR              << std::endl
#endif
	  << "Placchetta\t"     << nspt_pars.plaqn   << std::endl
#ifdef __WILLOOP_2x2__
	  << "Altri loops\t"    << nspt_pars.name2x2 << std::endl
#endif
	  << "Configurazione\t" << nspt_pars.confn   << std::endl
	  << "Damocle\t\t"      << nspt_pars.damon   << std::endl
	  << "Norma\t\t"        << nspt_pars.normn   << std::endl
	  << "Tr(U)\t\t"        << nspt_pars.trun    << std::endl << std::endl;

#ifdef __PARALLEL_OMP__

  cout << "\tThreads partitionning:" << endl;
	
  if ( (act_pars.sz[0]/ntt) < 3 ) { cout << "Errore: troppi threads nella direzione 0 rispetto alla taglia di reticolo." << endl; exit(-1); }
  if ( (act_pars.sz[1]/ntx) < 3 ) { cout << "Errore: troppi threads nella direzione 1 rispetto alla taglia di reticolo." << endl; exit(-1); }
  if ( (act_pars.sz[2]/nty) < 3 ) { cout << "Errore: troppi threads nella direzione 2 rispetto alla taglia di reticolo." << endl; exit(-1); }
  if ( (act_pars.sz[3]/ntz) < 3 ) { cout << "Errore: troppi threads nella direzione 3 rispetto alla taglia di reticolo." << endl; exit(-1); }
	
  for( int tid = 0; tid < NTHR; tid++){
    thr_pars[tid].xi = new int[dim];
    thr_pars[tid].xf = new int[dim];

    thr_pars[tid].xi[0] = ( act_pars.sz[0] / ntt ) * (  tid / (ntx*nty*ntz)       );
    thr_pars[tid].xi[1] = ( act_pars.sz[1] / ntx ) * ( (tid / (nty*ntz)   ) % ntx );
    thr_pars[tid].xi[2] = ( act_pars.sz[2] / nty ) * ( (tid / (ntz)       ) % nty );
    thr_pars[tid].xi[3] = ( act_pars.sz[3] / ntz ) * (  tid                 % ntz );
    
    thr_pars[tid].xf[0] = thr_pars[tid].xi[0] + act_pars.sz[0]/ntt ;
    thr_pars[tid].xf[1] = thr_pars[tid].xi[1] + act_pars.sz[1]/ntx ;
    thr_pars[tid].xf[2] = thr_pars[tid].xi[2] + act_pars.sz[2]/nty ;
    thr_pars[tid].xf[3] = thr_pars[tid].xi[3] + act_pars.sz[3]/ntz ;

    cout << "tid = " << tid << "\t"
	 << "\tseed = " << seed[tid] 
	 << std::endl
	 << thr_pars[tid].xi[0] << "\t"
	 << thr_pars[tid].xi[1] << "\t"
	 << thr_pars[tid].xi[2] << "\t"
	 << thr_pars[tid].xi[3] << "\n"
	 << thr_pars[tid].xf[0] << "\t"
	 << thr_pars[tid].xf[1] << "\t"
	 << thr_pars[tid].xf[2] << "\t"
	 << thr_pars[tid].xf[3] << "\n\n";
      
  }

#endif


#ifdef __K_MOM_ANALYSIS__
  char *kname = new char[100];
#ifndef __LANDAU_GAUGE_FIXING__
  sprintf(kname,"knorm_%s_nf0_%dx%dx%dx%d_a%2.2f_tg%3.3f_o%d.dat",
	  type_g, act_pars.sz[0], act_pars.sz[1], 
	  act_pars.sz[2],act_pars.sz[3],fabs(act_pars.alpha), act_pars.tau_g, 
	  PTORD);
#else
  sprintf(kname,"knorm_acc_%s_nf0_%dx%dx%dx%d_a%2.2f_tg%3.3f_o%d.dat",
	  type_g, act_pars.sz[0], act_pars.sz[1], 
	  act_pars.sz[2],act_pars.sz[3],fabs(act_pars.alpha), act_pars.tau_g, 
	  PTORD);
#endif
  FPkn = fopen(kname,"wb");
  delete [] kname;
#endif

#ifdef __TIMING__
  Time.init();
#endif

  return 0;
}





int QuenchedAllocate(ptGluon_fld& Umu){

  // Calcola il volume e la normalizzazione 1/V
  act_pars.iVol = ( act_pars.sz[0] * act_pars.sz[1] *
		    act_pars.sz[2] * act_pars.sz[3] );
  act_pars.rVol = 1.0/(double)act_pars.iVol;

  // Definisce sqrt(tau_g) e tau_g/(2Nc)
  act_pars.stau = sqrt(act_pars.tau_g); 
#ifdef _SQRT_BETA_EXPANSION_
  act_pars.tau_g *= -D6;
#elif defined _G_EXPANSION_
  act_pars.tau_g *= -1;
#endif

  // Genera il campo gluone freddo o da configurazione
  // secondo il flag nspt_pars.Init. Se legge da configurazione
  // legge anche la placchetta ad file e la confronta con una
  // misura fatta al volo

  if(!nspt_pars.Init){

    for(int i = 0; i < act_pars.iVol; i++){
      for(int mu = 0; mu < dim; mu++){
	Umu.W[i][mu].id();
      }
    }

  }
  else{
#ifndef APE_CONFIG    
    if( Umu.load(nspt_pars.confn, w1) ){
      std::cout << "Errore nella lettura della configurazione." << std::endl;
      return 1;
    }
#else
    Umu.load_ape();
#endif
    plaquette_measure(Umu, act_pars);
#ifndef APE_CONFIG
    if( plaquette_check(w, w1) ){
      std::cout << "La misura della placchetta non corrisponde." << std::endl;
      for (int i1=0; i1 <= PTORD; i1++){
	printf("%e %e\t",w[i1].re, w[i1].im );
	printf("%e %e\n",w1[i1].re,w1[i1].im);
      }
      //      return 1;
    }
#endif
    std::cout << "\nMisura della placchetta:" << std::endl;
    for (int i1=0; i1 <= PTORD; i1++){
      printf("%e %e\n",w[i1].re,w[i1].im);
    }

  }

  
  U = Umu.handle();

  // Data e ora di inizio run su logfile
  time_t tempo;
  time(&tempo);
  logfile << "Inizio run\t\t" << ctime(&tempo) << std::endl;

#ifdef __K_MOM_ANALYSIS__
  for(int i1 = 0; i1 <= allocORD;i1++)
    knorm[i1] = new double[Umu.Z->Size];
#endif

  return 0;
}



int AggiornaParametri(ptGluon_fld& Umu) {

#ifdef __WILLOOP_2x2__
  WL2x2(Umu);
#endif

#ifdef __MANY_LOOP__
  ComputeLoops(Umu);
#endif

  // Scrive log
  time_t tempo;
  time(&tempo);
  nInt += nspt_pars.Beat;
  logfile << "______________________________________"  << std::endl 
	  << std::endl;
  logfile << "Damocle Update\t\t" << ctime(&tempo)  << std::endl;
  logfile << "# Iterazioni\t\t"   << nInt           << std::endl;

  // Legge file damocle e aggiorna parametri
  FILE *df;
  df = fopen(nspt_pars.damon,"r");

  get_val(df, "SWEEP",       "%d" ,&(nspt_pars.Sweep) );
  get_val(df, "BEAT",        "%d" ,&(nspt_pars.Beat)  );
  get_val(df, "tau_g",       "%lf",&(act_pars.tau_g)  );
  get_val(df, "alpha",       "%lf",&(act_pars.alpha)  );
  get_val(df, "save_config", "%d" ,&(nspt_pars.Save)  );
  get_val(df, "kill_run",    "%d" ,&(nspt_pars.Kill)  );


  logfile << "Tau_g\t\t\t"        << act_pars.tau_g << std::endl;
  logfile << std::endl;
  logfile.flush();

  act_pars.stau = sqrt(act_pars.tau_g);
#ifdef _SQRT_BETA_EXPANSION_
  act_pars.tau_g *= -D6;
#elif defined _G_EXPANSION_
  act_pars.tau_g *= -1;
#endif

  // Se si salva il campo come checkpoint (save_config = 1)
  if(nspt_pars.Save == 1) {
    plaquette_measure(Umu, act_pars);
    Umu.save(nspt_pars.confn, w);
  }

  // Se si vuole terminare il run (kill_run = 1)
  if( nspt_pars.Kill == 1) { 
    nspt_pars.Sweep = 0;
  }

  return 0;
}


int NsptFinalize(ptGluon_fld& Umu, int t){

#ifdef __K_MOM_ANALYSIS__
  fclose(FPkn);
#endif

  time_t tempo; 
  time(&tempo);
  logfile << "__________________________________"  << std::endl
	  << "Termine simulazione\t\t"             << ctime(&tempo) << std::endl
	  << "# Iterazioni\t\t"                    <<  t << std::endl;

  plaquette_measure(Umu, act_pars);

  printf("Misura della placchetta:\n");
  for (int i1=0; i1 <= PTORD; i1++){
    printf("%e\t%e\n", w[i1].re, w[i1].im);
  }


#ifdef __WILLOOP_2x2__
  WL2x2(Umu);
  loop2x2.close();
#endif

#ifdef __MANY_LOOP__
  ComputeLoops(Umu);
  loop2x2.close();
#endif
  

  Umu.save(nspt_pars.confn,w);
  delete [] w;
  delete [] w1;
  delete [] w2;
  delete [] norm;
  delete [] norm1;

  return 0;
}



int main(int argc, char** argv){

  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);  // Enable fpe

  if( initialize(argc, argv) ){ exit(1); }


  // Alloca il reticolo e genera i momenti. A seconda del
  // flag ABC o PBC il reticolo sara' simmetrico o antisimmetrico
  // nella coordinata T
  latt LL(act_pars.sz);  
  LL.p_init();


  // Alloca il campo gluone (pt), fermione (pt) e un
  // campo fermionico ausiliario (non pt)
  ptGluon_fld Umu(&LL);

#ifdef __K_MOM_ANALYSIS__

#ifdef  __PARALLEL_OMP__
  fftw_init_threads();
  fftw_plan_with_nthreads(NTHR);
#endif

  planGluon = new fftw_plan[2];
  planGluon[0] = fftw_plan_many_dft(dim, LL.Sz, dim*(NC*NC*allocORD+1),
				     (fftw_complex *) Umu.W, 
				     NULL, dim*(NC*NC*allocORD+1), 1,
				     (fftw_complex *) Umu.W, 
				     NULL, dim*(NC*NC*allocORD+1), 1,
				     FFTW_FORWARD, FFTW_ESTIMATE);
  planGluon[1] = fftw_plan_many_dft(dim, LL.Sz, dim*(NC*NC*allocORD+1),
				     (fftw_complex *) Umu.W, 
				     NULL, dim*(NC*NC*allocORD+1), 1,
				     (fftw_complex *) Umu.W, 
				     NULL, dim*(NC*NC*allocORD+1), 1,
				     FFTW_BACKWARD, FFTW_ESTIMATE);
#endif

#ifdef __LANDAU_GAUGE_FIXING__

#ifdef  __PARALLEL_OMP__
  fftw_init_threads();
  fftw_plan_with_nthreads(NTHR);
#endif

  ptSU3_fld WFA(&LL);
  Wgauge = &WFA;

  planFA = new fftw_plan[2];
  planFA[0] = fftw_plan_many_dft(dim, LL.Sz, NC*NC*PTORD+1,
  				     (fftw_complex *) WFA.W, 
  				     NULL, NC*NC*allocORD+1, 1,
  				     (fftw_complex *) WFA.W, 
  				     NULL, NC*NC*allocORD+1, 1,
  				     FFTW_FORWARD, FFTW_ESTIMATE);
  planFA[1] = fftw_plan_many_dft(dim, LL.Sz, NC*NC*PTORD+1,
  				     (fftw_complex *) WFA.W, 
  				     NULL, NC*NC*allocORD+1, 1,
  				     (fftw_complex *) WFA.W, 
  				     NULL, NC*NC*allocORD+1, 1,
  				     FFTW_BACKWARD, FFTW_ESTIMATE);
#endif
  
  if( QuenchedAllocate(Umu) ) { exit(1); }

  // Se tauf_f > tau_g determina quante iterazioni gluoniche
  // tra due fermioniche
  int t1 = 0;

  while( t1 < nspt_pars.Sweep) {
    t1 += nspt_pars.Beat;

    NsptEvolve(Umu);
    
    // Legge i nuovi parametri da damocle
    if( AggiornaParametri(Umu) ){
      return 0;
    }
  }

  NsptFinalize(Umu, t1);
  return 0;
}
