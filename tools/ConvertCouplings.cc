#include"../QCDenvNODEpt.h"

#include<iostream>
#include<fstream>

extern int PTORD;

using namespace std;

int main(){

  PTORD = allocORD;

  int way;
  int *sz = new int[4];
  int *xx = new int[4];
  char *junk    = new char[100];
  char *confIn  = new char[100];
  char *confOut = new char[100];

  Cplx *plaq = new Cplx[PTORD];

  ifstream inFile;
  if (!inFile.is_open())
    {
      inFile.open("Conversion.txt");
      if (!inFile.is_open())
	{
	  cout << "Error reading configuration file." << endl;
	  return 1;
	}
    }

  inFile >> junk >> sz[0] >> sz[1] >> sz[2] >> sz[3];
  if( strcmp(junk,"Size") != 0) {
    cout << "Errore: atteso Size, trovato " << junk << endl;
    exit(0);
  }
  inFile >> junk >> confIn;
  if( strcmp(junk,"ConfIn") != 0) {
    cout << "Errore: atteso ConfIn, trovato " << junk << endl;
    exit(0);
  }
  inFile >> junk >> confOut;
  if( strcmp(junk,"ConfOut") != 0) {
    cout << "Errore: atteso ConfOut, trovato " << junk << endl;
    exit(0);
  }
  inFile >> junk >> way;
  if( strcmp(junk,"Way") != 0) {
    cout << "Errore: atteso Way, trovato " << junk << endl;
    exit(0);
  }

  cout << "Size = "
       << sz[0]     << "\t"
       << sz[1]     << "\t"
       << sz[2]     << "\t"
       << sz[3]     << endl
       << "ConfIn\t"  << confIn  << endl
       << "ConfOut\t" << confOut << endl << endl;
  cout << "Conversion:" << endl;
  if(way == 1)
    cout << "\t1/sqrt(beta) --> g" << endl;
  else
    cout << "\tg --> 1/sqrt(beta)" <<  endl;


  latt LL(sz);  

  ptGluon_fld Umu(&LL);
  Umu.load(confIn, plaq);
 
  double factor;
  if( way == 1)
    factor = 1/sqrt(2.0*NC);
  else
    if ( way == -1)
      factor = sqrt(2.0*NC);
  
  double *mult;
  mult = new double[PTORD];

  mult[0] = factor;
  if( way == 1)
    {
      cout << endl 
	   << "Order 0" << "\t"
	   << "1"
	   << "\t\t--> 1/sqrt(2Nc)^(0) = 1" << endl
	   << "Order 1" << "\t"
	   << mult[0] 
	   << "\t--> 1/sqrt(2Nc)^(1/2)" << endl;
      
      for( int ord = 1; ord < PTORD; ord++)
	{
	  mult[ord] = factor*mult[ord-1];
	  cout << "Order " << ord+1 << "\t"
	       << mult[ord] << "\t"
	       << "--> 1/sqrt(2Nc)^(" << ord+1
	       << "/2)" << endl;
	}
      cout << endl << endl;
    }
  else
    {
      cout << endl 
	   << "Order 0" << "\t"
	   << "1"
	   << "\t\t--> sqrt(2Nc)^(0) = 1" << endl
	   << "Order 1" << "\t"
	   << mult[0] 
	   << "\t\t--> sqrt(2Nc)^(1/2)" << endl;

      for( int ord = 1; ord < PTORD; ord++)
	{
	  mult[ord] = factor*mult[ord-1];
	  cout << "Order " << ord+1 << "\t"
	       << mult[ord] << "\t"
	       << "\t--> sqrt(2Nc)^(" << ord+1
	       << "/2)" << endl;
	}
      cout << endl << endl;
      
    }
  
  
  for (int i= 0;i < LL.Size; i++){
   
    for( int mu = 0; mu < 4; mu++){

      for( int ord = 0; ord < PTORD; ord++){
	

  	Umu.W[i].U[mu].ptU[ord] = mult[ord] * Umu.W[i].U[mu].ptU[ord];


      }
 
    }


  }





  
  for(int i1 = 0; i1 < PTORD+1; i1++){
    plaq[i1] = 0.0;
  }					
  
  ptSU3 W1, W2;
  W1.zero();				  
  
  for(int i = 0; i < LL.Size;i++){

    for(int mu = 0; mu < dim; mu++){
      W2.zero();
      
      for(int nu = 0; nu < dim; nu++){
	
	if(nu != mu ){
	  
	  W2 += Umu.staple(i, mu, nu);
	  
	}
	
      }
      
      W1 += Umu.W[i].U[mu]*W2;
      
    }	
  
  }
  
  
  W1.Tr(plaq);
  
  cout << "Plaquette" << endl;
  for(int i1 = 0; i1 <= PTORD; i1++){  
    plaq[i1] /= (LL.Size*72);
    plaq[i1].prout();
    cout << endl;
  }					




  Umu.save(confOut, plaq);


  delete [] xx;
  delete [] sz;
  
  inFile.close();

  return 0;

}
