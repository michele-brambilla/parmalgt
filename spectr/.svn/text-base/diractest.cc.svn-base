#include <unistd.h>
#include<fstream>

#include <MassC.h>
#include <dirac_spec.h>
#include <dirac_par.h>

int mo; // Perturbative order required for corrections. 
int fs; // First degeneration subspace to be considered.
int ls; // Last degeneration subspace to be considered (0 means last one).
int vb; // Absolute value => How much output. Positive => Human formatted. Negative => Raw formatted.
int rw;
int lg;

int *size;
char *leadconf;
char *umuconf; 
char *outfile;


using namespace std;

int main(int argc, char **argv) {

  //  cout.precision(5);

  size     = new int[dim];
  leadconf = new char[100];
  umuconf  = new char[100];
  outfile  = new char[100];
  int sz, m, ns=0;
  eigv *ev;
  dgsp *ds;
  int* xx = new int[dim];

  // Legge i parametri
  set_params(argc, argv);

  // Alloca i campi e lo spazio di lavoro
  latt Z(size);
  Z.p_init();

  sz = Z.Size * dim * NC;
  ev = (eigv *) malloc(sz * sizeof(eigv));


//   // Genera l'autovalore libero ?? posso evitare dws.dd ??
//   for (int i = 0; i < z.Size; i++){
//     for (int mu = 0; mu < dim; mu++){
//       for (int a = 0; a < NC; a++) {
// 	m = i*NC*dim + mu*NC + a;
// 	dws.p->ket(i, mu, a);
// 	dws.dd(umu, 0);
// 	ev[m].id.i = i;
// 	ev[m].id.mu = mu;
// 	ev[m].id.a = a;
// 	ev[m].mo = 0;
// 	ev[m].pt[0] = ((dws.p)->bra(i, mu, a)).re;
//       }
//     }
//   }

  // Genera l'autovalore libero
  for (int i = 0; i < Z.Size; i++){
    Z.get(i,xx);
    double eval = MBARE + .5*(Z.p2hat[0][xx[0]] + Z.p2hat[1][xx[1]] +
			      Z.p2hat[2][xx[2]] + Z.p2hat[3][xx[3]]);
    eval = eval*eval + ( Z.pbar[0][xx[0]]*Z.pbar[0][xx[0]] +
			 Z.pbar[1][xx[1]]*Z.pbar[1][xx[1]] +
			 Z.pbar[2][xx[2]]*Z.pbar[2][xx[2]] +
			 Z.pbar[3][xx[3]]*Z.pbar[3][xx[3]] );

    for (int mu = 0; mu < dim; mu++){
      for (int a = 0; a < NC; a++) {
	m = i*NC*dim + mu*NC + a;
	ev[m].id.i = i;
	ev[m].id.mu = mu;
	ev[m].id.a = a;
	ev[m].mo = 0;
	ev[m].pt[0] = eval;
      }
    }
  }



  ptGluon_fld umu(&Z);
  umu.load(umuconf);

  diracws dws(&Z);
    
  conta_sottospazi(ev, sz, fs, ls, ns, 0);

  ds = new dgsp[ns];

  degen_subspc(ev, ds, sz, fs, ls, ns, 0);

  if (fs > ns) fs = ns-1;
  if ((ls > ns) || ls == 0) ls = ns;

  Cplx x;

  for(int r = fs; r < ls; r++){
    cout << "\tr = " << r << endl;
    cout << "\tr = " << ds[r].sz << endl;
    
    eigsysws ews(&Z, ds[r].sz);    
    
    for(int i = 0; i < ds[r].sz; i++){
      Z.get(ev[i+ds[r].bs].id.i,xx);
//       cout << xx[0] << "\t"
// 	   << xx[1] << "\t"
// 	   << xx[2] << "\t"
// 	   << xx[3] << endl;
      
      dws.p->zeros();
      dws.p->psi[Z.get(ev[i+ds[r].bs].id.i)].psi[ev[i+ds[r].bs].id.mu].whr[ev[i+ds[r].bs].id.a].re = 1;
      dws.dd(umu, 1);

      for(int j = 0; j < ds[r].sz; j++){
	x = dws.p->psi[Z.get(ev[j+ds[r].bs].id.i)].psi[ev[j+ds[r].bs].id.mu].whr[ev[j+ds[r].bs].id.a];
	//	cout << i << "\t" << j << endl;
	gsl_matrix_complex_set(ews.matrix, j, i, gsl_complex_rect(x.re,x.im));
      }
    }
    //return 0;


    ews.prout(ews.matrix);
    return 0;    
    ews.diagonalize();
    gsl_eigen_genhermv_sort(ews.eigval, ews.eigvec, GSL_EIGEN_SORT_VAL_ASC);
    
    for(int p = ds[r].bs; p < ds[r].bs+ds[r].sz; p++){
      ev[p].mo = 1;
      ev[p].pt[1] = ews.eigval_get(p-ds[r].bs);
      cout << ev[p].pt[1] << endl;
    }



    dws.p->zeros();
    for(int i = ds[r].bs; i < ds[r].bs+ds[r].sz; i++){
      for(int j = ds[r].bs; j < ds[r].bs+ds[r].sz; j++){
	dws.p->psi[Z.get(ev[j].id.i)].psi[ev[j].id.mu].whr[ev[j].id.a] = 
	  ews.eigvec_get(j,i);
      }
    }
    SpinColor_fld pp(&Z);
    pp = *dws.p;

    dws.dd(umu,1);
    //    cout << "----------------"<< endl;
//     for(int i = ds[r].bs; i < ds[r].bs+ds[r].sz; i++){
//     cout << ev[j].id.i << "\t"
// 	 << ev[j].id.mu << "\t"
// 	 << ev[j].id.a << "\t"
// 	 << ev[j].pt[1] << "\t";
    
//     for(int j = ds[r].bs; j < ds[r].bs+ds[r].sz; j++){
//       dws.p->psi[Z.get(ev[j].id.i)].psi[ev[j].id.mu].whr[ev[j].id.a].prout();
//       cout << "\t";
//       pp.psi[Z.get(ev[j].id.i)].psi[ev[j].id.mu].whr[ev[j].id.a].prout();
//       cout << "\t";
//       (dws.p->psi[Z.get(ev[j].id.i)].psi[ev[j].id.mu].whr[ev[j].id.a] - 
//        ev[j].pt[1]*pp.psi[Z.get(ev[j].id.i)].psi[ev[j].id.mu].whr[ev[j].id.a]).prout();
//       cout << endl;
//     }
//     }


    conta_sottospazi(ev,ds[r].sz, ds[r].bs, ds[r].bs+ds[r].sz, m, 1);

    ds[r].ds = new dgsp[m];
    
    cout << "r = " << r << "\tdimensione = " << ds[r].sz
	 << "\tnumero sottospazi= " << m << endl;


    if( (ds[r].ns=m) != ds[r].sz){
      cout << "Ulteriore degenerazione" << endl;
      degen_subspc(ev, ds[r].ds, ds[r].sz, ds[r].bs, 
		   ds[r].bs+ds[r].sz, m, 1);

//       for(int s = 0; s < m; s++){
// 	for(int p = 0; p < ds[r].ds[s].sz; p++){
// 	  cout << "\ts = "  << s 
// 	       << "\ti = "  << ev[ds[r].ds[s].bs+p].id.i
// 	       << "\tmu = " << ev[ds[r].ds[s].bs+p].id.mu
// 	       << "\ta = "  << ev[ds[r].ds[s].bs+p].id.a
// 	     << endl;
// 	}
// 	cout << "\tr = "<< r << "\ts = " << s 
// 	     << "\tbs = " << ds[r].ds[s].bs
// 	     << "\tsz = " << ds[r].ds[s].sz
// 	     << "\t" << ds[r].ev
// 	     << "\t" << ds[r].ds[s].ev
// 	     << endl;
//       }
      
      for(int s = 0; s < m; s++){
	eigsysws ews1(&Z, ds[r].ds[s].sz);	
	
	for(int i = ds[r].ds[s].bs; i < ds[r].ds[s].bs+ds[r].ds[s].sz; i++){

	  dws.p->zeros();

	  for(int j = ds[r].ds[s].bs; j < ds[r].bs+ds[r].sz; j++){
	    dws.p->psi[Z.get(ev[j].id.i)].psi[ev[j].id.mu].whr[ev[j].id.a] =
	      ews.eigvec_get(i-ds[r].ds[s].bs, j-ds[r].bs );
	  }

	  return 0;
	  dws.dd(umu, 2);

	  for(int j = ds[r].ds[s].bs; j < ds[r].ds[s].bs+ds[r].ds[s].sz; j++){
	    
 	    x = dws.p->psi[ev[j].id.i].psi[ev[j].id.mu].whr[ev[j].id.a];
 	    
	    gsl_matrix_complex_set(ews1.matrix, i-ds[r].ds[s].bs, 
				   j-ds[r].ds[s].bs, gsl_complex_rect(x.re,x.im));
	  }
 	
	}
	//	ews1.prout(ews1.matrix);
	//ews.diagonalize();
	
	for(int i = ds[r].ds[s].bs; i < ds[r].ds[s].bs+ds[r].ds[s].sz; i++){
	  ev[i].mo = 2;
	  ev[i].pt[1] = gsl_vector_get(ews.eigval, i-ds[r].ds[s].bs);
	}
	
	
	}
      
      
      //       for(int s = ds[r].ds[s].bs; s < ds[r].ds[s].bs+ds[r].ds[s].sz; s++){
// 	cout << ds[r].ns << "\t" << ev[s].pt[0] << "\t"
// 	     << ev[s].pt[1] << "\t" << ev[s].id.i << "\n";

	
// 	for(int i = ds[r].ds[s].bs; i < ds[r].ds[s].bs+ds[r].ds[s].sz; i++){

// 	  dws.p->psi[ev[i].id.i].psi[ev[i].id.mu].whr[ev[i].id.a].re = 1;
// 	  dws.dd(umu, 1);
// 	}
	
//  	eigenproblem(ews1, ds[r].ds[s], dws, ev, umu, 2);
// 	conta_sottospazi(ev, ds[r].ds[s].sz, ds[r].ds[s].bs, 
// 			 ds[r].ds[s].bs+ds[r].ds[s].sz, m, 2);
// 	ds[r].ds[s].ds = new dgsp[m];
	
// 	cout << "\t\tr= " << r 
// 	     << "\ts=" << s 
// 	     << "\tdimensione = " << ds[r].ds[s].sz
// 	     << "\tnumero sottospazi= " << (ds[r].ds[s].ns=m) << endl;
//       }

//       cout << "uscito" << endl;

/*    else{
      cout << "Degenerazione risolta - ordine " << ev[ds[r].bs].mo << endl;
      for(int i = 0; i < ds[r].sz; i++){
	ds[r].ds[i].bs = ds[r].bs;
	ds[r].ds[i].sz = ds[r].sz;
	ds[r].ds[i].ns = ds[r].ns;
      }
      

      non_degen_solver(dws, ev, ews, ds[r], umu, r, mo);

      for(int i = ds[r].bs; i < ds[r].bs+ds[r].sz; i++){
	cout << "\t" << ev[i].id.i;
	for(int k = 0; k < mo; k++){
	  cout << "\t" << ev[i].pt[k];
	}
	cout << endl;
      } */


    }
  
  
  }
  
  
  
  
  free(ev);
  
  return 0;
}

