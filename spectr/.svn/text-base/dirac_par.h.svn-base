#include<iostream>

extern int mo;
extern int fs;
extern int ls;
extern int vb;
extern int rw;
extern int lg;

extern int *size;
extern char *leadconf;
extern char *umuconf; 
extern char *outfile;



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




void set_params(int argc, char **argv){
  
  // Valori di default
  mo = 2;
  fs = 0;
  ls = 0;
  vb = 2;
  rw = 0;
  lg = 0;
  
  // Legge valori da file
  FILE *fp;
  if( (fp = fopen("Dirac.cfg", "r")) == NULL) { exit(1); }
  fscanf(fp, "size     %d %d %d %d\n", size, size+1, size+2, size+3);
  get_val(fp, "maxord"   ,"%d", &mo);
  get_val(fp, "firstsp"  ,"%d", &fs);
  get_val(fp, "lastsp"   ,"%d", &ls);
  get_val(fp, "verbose"  ,"%d", &vb);
  get_val(fp, "raw"      ,"%d", &rw);
  get_val(fp, "leadgen"  ,"%d", &lg);
  get_val(fp, "leadconf" ,"%s", leadconf);
  get_val(fp, "umuconf"  ,"%s", umuconf);
  get_val(fp, "outfile"  ,"%s", outfile);
  fclose(fp);

  // I comandi da terminale dominano
  int opt;
  extern char *optarg;
  while ((opt = getopt(argc, argv, "f:l:o:r:v:")) != -1) {
    switch (opt) {
    case 'f':
      fs = atoi(optarg);
      break;
    case 'l':
      ls = atoi(optarg);
      break;
    case 'o':
      mo = atoi(optarg);
      break;
    case 'r':
      rw = atoi(optarg);
      break;
    case 'v':
      vb = atoi(optarg);
      break;
    case '?':
      fprintf(stderr, "Unrecognized option.\n");
      break;
    }
  }

  printf( "Taglia   = %d %d %d %d\n", size[0], size[1], size[2], size[3]);
  printf( "maxord   = %d\n", mo);
  printf( "firstsp  = %d\n", fs);
  printf( "lastsp   = %d\n", ls);
  printf( "verbose  = %d\n", vb);
  printf( "raw      = %d\n", rw);
  printf( "leadgen  = %d\n", lg);
  printf( "leadconf = %s\n", leadconf);
  printf( "umuconf  = %s\n", umuconf);
  printf( "outfile  = %s\n", outfile);


}
