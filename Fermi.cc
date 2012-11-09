#include <LocalField.hpp>
#include <Background.h>
#include <Geometry.hpp>
#include <algorithm>
#include <newQCDpt.h>
#include <map>
#include <iostream>
#include <fstream>
#include <Kernels.hpp>
#include <uparam.hpp>
#include <stdlib.h>
#include <IO.hpp>
#ifdef _OPENMP
#include <omp.h>
#else
#include <time.h>
#endif
#include <signal.h>

#ifdef USE_MPI
#include <mpi.h>
#include <sstream>
#endif

bool soft_kill = false;
int got_signal = 100;
void kill_handler(int s){
       soft_kill = true;
       got_signal = s;
       std::cout << "INITIATE KILL SEQUENCE\n"; 
}

// space-time dimensions
const int DIM = 4;
// perturbative order
const int ORD = 4;

// max number of iterations
const int COUNT_MAX = 100;


////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
//
// PARAMETERS to be read later

// lattice size
int L;
int T;
// s parameter for staggered
int s;
// total number of gauge updates
int NRUN;
// frequency of measurements
int MEAS_FREQ;
// testing gauge fixing option -- DO NOT TOUCH!
const int GF_MODE = 1;
// integration step and gauge fixing parameter
double taug;
double alpha;

// some short-hands
typedef bgf::ScalarBgf Bgf_t; // background field
// typedef bgf::AbelianBgf Bgf_t; // background field
typedef BGptSU3<Bgf_t, ORD> ptSU3; // group variables
typedef ptt::PtMatrix<ORD> ptsu3; // algebra variables
typedef BGptGluon<Bgf_t, ORD, DIM> ptGluon; // gluon
typedef pt::Point<DIM> Point;
typedef pt::Direction<DIM> Direction;

// shorthand for gluon field
typedef fields::LocalField<ptGluon, DIM> GluonField;
typedef GluonField::neighbors_t nt;

// shorthand for fermion field
typedef SpinColor<4> Fermion;
typedef fields::LocalField< Fermion , DIM> ScalarFermionField;
typedef std::vector<ScalarFermionField> FermionField;


//
// Make aliases for the Kernels ...
//

typedef kernels::WilsonPTKernel<GluonField> WilsonKernel;
typedef kernels::WilsonTreeLevelKernel<GluonField> WilsonTLKernel;


// ... to set the background field ...
typedef kernels::SetBgfKernel<GluonField> SetBgfKernel;



// timing

struct Timer {
  double t, tmp;
  static double t_tot;
  Timer () : t(0.0) { };
  void start() { 
#ifdef _OPENMP
    tmp = omp_get_wtime(); 
#else
    tmp = clock();
#endif

  }
  void stop() { 
#ifdef _OPENMP
    double elapsed = omp_get_wtime() - tmp; 
#else
    double elapsed = ((double)(clock() - tmp))/CLOCKS_PER_SEC;
#endif
    t += elapsed;
    t_tot += elapsed;
  }
};

double Timer::t_tot = 0.0;

// helper function to feed int, double etc. to the parameters
template <typename T> 
std::string to_string(const T& x){
  std::stringstream sts;
  sts << x;
  return sts.str();
}



// BCGstab
template<int DIM>
void inverter( GluonField& U, ScalarFermionField& x, double& m )
{

  const double inv_prec = 1e-5;
  Cplx alpha, beta, omega;
  Cplx nr, nr1;
  Cplx norm_r, norm_r1;

  //  b.randomize();  // ok, here preconditioning or something else...

  ScalarFermionField b (x), r0(x), Ax(x);
  WilsonTLKernel apply_from_x(U, x, m );
  
  Ax.apply_everywhere(apply_from_x);
  r0 = b - Ax;

  ScalarFermionField p(x), p0(r0), r(r0), r0star(r), s(x);
  WilsonTLKernel apply_from_p(U, p, m );
  WilsonTLKernel apply_from_s(U, s, m );

  ScalarFermionField Ap(x), As(x);

  int count = 0;

  nr1 = (r * r0star);
  while( count < COUNT_MAX )
    {
  
      ++count;
      Ap.apply_everywhere(apply_from_p);

      nr = nr1;

      Ap.apply_everywhere(apply_from_p);

      alpha = nr / ( Ap * r0star );
      s     = r - Ap * alpha;
      
      As.apply_everywhere(apply_from_s);

      omega = ( As * s ) / ( As * As );
      x     = x + ( (p * alpha) + (s * omega) );
      r     = s - As * omega;
      nr1   = (r * r0star);
      beta  = (nr1 / nr) * (alpha / omega);
      p     = r + (p - Ap * omega) * beta;
      
      std::cout << count << "\t" << r*r << "\n";
    }


}






void invertQ( ScalarFermionField& src, GluonField& U, FermionField& dest ) {
  
  std::cout << "Fermion Operator Inversion" << std::endl;
  
  double m = 0.0;
  array_t<double, ORD>::Type mass;
  
  // WilsonTLKernel wtlk(U, src, mass[0] );
  
  mass[0] = 1;

  // dest[0].apply_everywhere(wtlk);  
  // WilsonKernel wk(U, dest, mass );
  // for( wk.ord() = 1; wk.ord() < ORD; ++wk ) {    
  //   dest[wk.ord()].apply_everywhere(wk);
  // }
  
  inverter<DIM>( U, dest[0], mass[0] );

}





int main(int argc, char *argv[]) {

  signal(SIGUSR1, kill_handler);
  signal(SIGUSR2, kill_handler);
  signal(SIGXCPU, kill_handler);
  int rank;

  ////////////////////////////////////////////////////////////////////
  // read the parameters
  uparam::Param p;
  p.read("input");
  L = atof(p["L"].c_str());
  s = atoi(p["s"].c_str());
  alpha = atof(p["alpha"].c_str());
  taug = atof(p["taug"].c_str());
  NRUN = atoi(p["NRUN"].c_str());
  MEAS_FREQ = atoi(p["MEAS_FREQ"].c_str());
  T = L-s;
  // also write the number of space-time dimensions
  // and perturbative order to the parameters, to
  // make sure they are written in the .info file 
  // for the configurations stored on disk
  p["NDIM"] = to_string(DIM);
  p["ORD"] = to_string(ORD);
  ////////////////////////////////////////////////////////////////////
  //
  // timing stuff
  typedef std::map<std::string, Timer> tmap;
  tmap timings;
  ////////////////////////////////////////////////////////////////////
  //
  // random number generators
  srand(atoi(p["seed"].c_str()));
  ////////////////////////////////////////////////////////////////////
  //
  // lattice setup
  // generate an array for to store the lattice extents
  geometry::Geometry<DIM>::extents_t e;
  // we want a L = 4 lattice
  std::fill(e.begin(), e.end(), L);
  // for SF boundary: set the time extend to T + 1
  e[0] = T + 1;
  // we will have just one field
  GluonField U(e, 1, 0, nt());



  ScalarFermionField Xi  (e, 1, 0, nt());
  FermionField Psi;

  int n = 0;
  for ( ScalarFermionField::iterator it = Xi.begin(); it != Xi.end(); ++it, ++n )
    for( Direction mu(0); mu.is_good(); ++mu) 
      (*it)[mu].whr[0] = Cplx(n,1);

  for( int i = 0; i < ORD; ++i)
    Psi.push_back(ScalarFermionField(e, 1, 0, nt()));


  invertQ( Xi, U, Psi );

  return 0;
}
