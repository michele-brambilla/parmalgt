#include<dirac_spec.h>
//#include<fstream>

using namespace std;

int main(){
  int *sz, ptord, seed, rep;
  char conf[100];
  sz = new int(dim);
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
  SpinColor_fld Xi(&LL);
  ptSpinColor_fld Pmu(&LL);

  LL.p_init();

  Cplx **matr, res[ptord-2];

  matr = new Cplx*[ptord-2]; 
  for(int i = 0; i < ptord-2; i++)
    matr[i] = new Cplx[LL.Size*NC*dim];

  Umu.load(conf);
  rand_init(seed);

/*   ofstream(fout); */
/*   fout.open("test_dd.dat"); */

  for(int rp = 0; rp < rep; rp++){
    
    Xi.gauss();
    *(Pmu.scfld) = Xi;
    
    Pmu.fillPT(ptord-2, Umu);
    
    for(int ord = 0; ord < ptord-2; ord++){
      for(int i = 0; i < LL.Size; i++)
	for(int mu = 0; mu < dim; mu++)
	  for(int a = 0; a < NC; a++){
	    res[ord] += (matr[ord][a+NC*mu+dim*i]
			 = (~Xi.psi[i].psi[mu].whr[a]*Pmu.psi[i].psi[mu].ptCV[ord].whr[a]) );
/* 	    fout << matr[ord][a+NC*mu+dim*i].re << "\t" */
/* 		 << matr[ord][a+NC*mu+dim*i].im << "\t" */
/* 		 << endl; */
	  }
    }
    
    for(int i = 0; i < ptord-2; i++){
      cout << res[i].re/(double)(LL.Size*dim*NC) << "\t" 
	   << res[i].im/(double)(LL.Size*dim*NC) << "\t"
	   << endl;
      res[i] = 0;
    }
  
  }

  //    for(int i = 0; i < LL.Size*NC*dim; i++){
  //      matr[0][i].prout();
  //      cout << endl;
  //    }
  

  //  fout.close();
  delete [] matr;
  return(0);
}



