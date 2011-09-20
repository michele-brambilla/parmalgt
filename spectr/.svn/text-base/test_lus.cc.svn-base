#include<dirac_spec.h>
#include<fftw3.h>

using namespace std;


int main(){

  int *sz, ptord, seed, rep;
  sz = new int(dim);
  char conf[100];
  FILE *fin;

  fin = fopen("TestSpectr.cfg","r");
  fscanf(fin, "taglia %d %d %d %d \n", &sz[0], &sz[1], &sz[2], &sz[3]);
  fscanf(fin, "ordine %d \n",          &ptord);
  fscanf(fin, "seed %d \n",            &seed);
  fscanf(fin, "conf %s \n",            &conf);
  fscanf(fin, "ripetiz %d \n",         &rep);
  fclose(fin);


  latt LL(sz);
  ptGluon_fld Umu(&LL);
  ptSpinColor_fld Pmu(&LL);
  diracws dws(&LL);
  LL.p_init();

  Cplx  res[ptord+1], norm;

  Umu.load(conf);

  rand_init(seed);

  for(int rp = 0; rp < rep; rp++){
    
    Pmu.scfld->gauss();
    for(int ord = 0; ord <= ptord; ord++){
      *(dws.tq) = *(Pmu.scfld);
      
      dws.fft(1);
      dws.dd(Umu,ord);
      dws.fft(0);
      
      for(int i = 0; i < LL.Size; i++)
	for(int mu = 0; mu < dim; mu++)  
	  Pmu.psi[i].psi[mu].ptCV[ord] = dws.q->psi[i].psi[mu];
    }
    
    for(int ord = 0; ord <= ptord; ord++){
      for(int i = 0; i < LL.Size; i++)
	for(int mu = 0; mu < dim; mu++)
	  for(int a = 0; a < NC; a++){
	    res[ord] += (~(Pmu.scfld->psi[i].psi[mu].whr[a])*Pmu.psi[i].psi[mu].ptCV[ord].whr[a]) ;
	}
    }
    
    for(int ord = 0; ord <= ptord; ord++){
      cout << (res[ord]/(double)(LL.Size*dim*NC)).re << "\t" 
	   << (res[ord]/(double)(LL.Size*dim*NC)).im
	   << endl;
      res[ord] = 0;
    }
    
  }

  return(0);
}


