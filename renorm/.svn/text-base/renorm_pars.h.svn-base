#ifndef _REN_PAR_H_
#define _REN_PAR_H_

#include<iostream>
#include<unistd.h>

//typedef struct{
//  int iVol;
//  int ptord;
//  double alpha;
//  int *sz;
//  char *conf;
//  char *out; 
//  char *prfile;
//  char *prfileT;
//} renorm_params_t;


using namespace std;

extern renorm_params_t params;

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





int get_params(int argc, char **argv){

  FILE *fin;
  if( (fin = fopen("Ren.cfg","r")) == NULL){
    return 1;
  }
  else{

    fscanf(fin, "taglia %d %d %d %d \n", &params.sz[0], &params.sz[1], &params.sz[2], &params.sz[3]);

    get_val(fin, "ordine", "%d",  &(params.ptord)  );
    get_val(fin, "conf"  , "%s",   (params.conf)   );
    get_val(fin, "out"   , "%s",   (params.out)    );
    get_val(fin, "prop"  , "%s",   (params.prfile) );
    get_val(fin, "alpha" , "%lf", &(params.alpha)  );
    //    get_val(fin, "nthr"  , "%d",  &(nthr)          );

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


#ifdef __TRANSP_PROPAG__
    strncpy(params.prfileT,params.prfile,strlen(params.prfile)-4);
    strcat(params.prfileT,".tra");
#endif

    printf("%s\n",params.conf);

    cout << "Parametri:\n" 
	 << "taglia " <<params.sz[0]<<" "<< params.sz[1] <<" "<< params.sz[2] <<" "<< params.sz[3] <<"\n"
	 << "ordine "         << params.ptord   << "\n"
	 << "configurazione " << params.conf    << "\n"
	 << "output "         << params.out     << "\n"
	 << "propagatore "    << params.prfile  << "\n"
#ifdef __TRANSP_PROPAG__
	 << "propagatoreTra " << params.prfileT << "\n"
#endif
	 << "alpha "          << params.alpha   << "\n"
      //	 << "nthr "           << nthr           << "\n"
	 /* << "R "              << params.R       << "\n\n"; */
	 << endl;
  }
  
  return 0;
}


#endif
