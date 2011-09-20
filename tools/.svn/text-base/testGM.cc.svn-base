#include<GammaStuff.h>
#include<MyQCD.h>

using namespace std;


extern Cplx gmuval[15][4];
extern int gmuind[15][4];







int main(int argc, char **argv){

  MyRand rand(12345);
  SpinColor p, gmp;

  for(int mu = 0; mu < dim; mu++)
    {
      for(int a = 0; a < NC; a++)
	{
	  //	  p.psi[mu].whr[a] = Cplx(rand.Rand(),rand.Rand());
	  p.psi[mu].whr[a] = Cplx(mu+1,a+1);
	}
    }

  //S
  for(int cid = 0; cid < 4; cid++)
    {
      cout << "cid = "     << cid << endl;
      
      //p.gmleft(cid).prout();
      //cout << "--------------------" << endl;
      for(int nu = 0; nu < dim; nu++)
	gmp.psi[nu] = gmuval[cid][nu]*(p.psi[gmuind[cid][nu]]);
      
      (gmp-p.gmleft(cid)).prout();
    }


  //P
  cout << "cid = "     << 5 << endl;
  
  //p.gmleft(cid).prout();
  //cout << "--------------------" << endl;
  for(int nu = 0; nu < dim; nu++)
    gmp.psi[nu] = gmuval[4][nu]*(p.psi[gmuind[4][nu]]);
      
  (gmp-p.gmleft(5)).prout();
  

  //A
  for(int cid = 0; cid < 4; cid++)
    {
      cout << "cid = "     << cid+5 << endl;
      
      //p.gmleft(cid).prout();
      //cout << "--------------------" << endl;
      for(int nu = 0; nu < dim; nu++)
	gmp.psi[nu] = gmuval[5+cid][nu]*(p.psi[gmuind[5+cid][nu]]);
      
      (gmp-p.gmleft(5).gmleft(cid)).prout();
    }


  //T
  int id = 9;
  for(int cid = 0; cid < 4; cid++)
    {
      for(int fid = cid+1; fid < 4; fid++)
	{
	  cout << "id = "     << cid << fid << endl;
      
	  for(int nu = 0; nu < dim; nu++)
	    gmp.psi[nu] = gmuval[id][nu]*(p.psi[gmuind[id][nu]]);
	  
	  (gmp-p.gmleft(fid).gmleft(cid)).prout();
	  id++;
	}
    }





}

