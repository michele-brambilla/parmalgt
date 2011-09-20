#include"dirac_spec.h"
#include<fftw3.h>

using namespace std;

int main(){
  int *sz, ptord, seed;
  char conf[100];
  sz = new int(dim);
  FILE *fin;

  fin = fopen("testconf.txt","r");

  fscanf(fin, "taglia %d %d %d %d \n", &sz[0], &sz[1], &sz[2], &sz[3]);
  fscanf(fin, "ordine %d \n",          &ptord);
  fscanf(fin, "seed %d \n",            &seed);
  fscanf(fin, "conf %s \n",            &conf);
  fclose(fin);

  latt LL(sz);
  ptGluon_fld Umu(&LL);
  ptSpinColor_fld Pmu(&LL, &Umu);
  Cplx **matr, res[ptord-2];

  matr = new Cplx*[ptord-2]; 
  for(int i = 0; i < ptord-2; i++)
    matr[i] = new Cplx[LL.Size*NC*dim];

  Umu.load(conf);

  rand_init(seed);
  Pmu.xi->gauss();

  Pmu.fillPT(ptord-2);

  for(int ord = 0; ord < ptord-2; ord++){
    for(int i = 0; i < LL.Size; i++)
      for(int mu = 0; mu < dim; mu++)
	for(int a = 0; a < NC; a++){
	  res[ord] += (matr[ord][a+NC*mu+dim*i]
		       = (~Pmu.xi->psi[i].psi[mu].whr[a]*Pmu.psi[i].psi[mu].ptCV[ord].whr[a]) );
	  matr[ord][a+NC*mu+dim*i].prout();
	  cout << endl;
	}
  }

  cout << endl;
  res[0].prout();
  cout << endl;
  res[1].prout();
  cout << endl;
//    for(int i = 0; i < LL.Size*NC*dim; i++){
//      matr[0][i].prout();
//      cout << endl;
//    }


  
  return(0);
}



