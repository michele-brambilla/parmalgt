/*
testAccuracy: testAccuracy.o MyQCD.o QCDpt.o MyMath.o MyRand.o
	$(CC) $(MyO)  -o testAccuracy.exe testAccuracy.o MyQCD.o QCDpt.o MyMath.o MyRand.o -L$(FFTWdir)lib -lfftw3 $(LibThr) -I$(FFTWdir)include

testAccuracy.o: testAccuracy.cc lattice.h QCDenvNODEpt.h QCDenvNODE.h
	$(CC) $(MyO) -c testAccuracy.cc -I. -I$(FFTWdir)include
*/

#include <iostream>
#include "../QCDenvNODEpt.h"

using namespace std;

int main(int argc, char **argv)
{

  int *sz = new int[dim];

  sz[0] = atoi(argv[1]);
  sz[1] = atoi(argv[1]);
  sz[2] = atoi(argv[1]);
  sz[3] = atoi(argv[1]);

  latt LL(sz);
  ptGluon_fld Umu(&LL);
  ptSU3 tmp;


  Umu.load(argv[2]);

  int i= 1;

  while(i!=0)
    {
      cin >> i;
	
      Umu.W[i].U[0].prout();
      cout << endl;
      tmp = log(Umu.W[i].U[0]);
      (Umu.W[i].U[0] - exp(tmp)).prout();

      cout << endl;
      cout << endl;

    }

  return 0;


}
