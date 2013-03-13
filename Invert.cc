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
#include <newMyQCD.h>
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
const int COUNT_MAX = 1000;


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

enum bdy{bulk=0,lower=1,upper=2};
template <class T>
struct FermionBulkKernels
{
  typedef kernels::WilsonTreeLevel5Kernel<T,bulk> WilsonTLKernel;
  typedef kernels::Wilson5Kernel<T,bulk> WilsonKernel;
};

template <class T>
struct FermionLowerKernels
{
  typedef kernels::WilsonTreeLevel5Kernel<T,lower> WilsonTLKernel;
  typedef kernels::Wilson5Kernel<T,lower> WilsonKernel;
};

template <class T>
struct FermionUpperKernels
{
  typedef kernels::WilsonTreeLevel5Kernel<T,upper> WilsonTLKernel;
  typedef kernels::Wilson5Kernel<T,upper> WilsonKernel;
};

typedef FermionBulkKernels<GluonField>  Bulk;
typedef FermionLowerKernels<GluonField> Lower;
typedef FermionUpperKernels<GluonField> Upper;


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
template<int DIM, int count_max>
int BCGstab0(  GluonField& U,  ScalarFermionField& src, ScalarFermionField& x,  double& m, const double& tol )
{

  x = src;
  Cplx alpha, beta, omega;
  Cplx nr, nr1;
  Cplx norm_r, norm_r1;

  //  b.randomize();  // ok, here preconditioning or something else...

  ScalarFermionField b (x), r0(x), Ax(x);
  //  WilsonTLKernel apply_from_x(U, x, m );
  Bulk::WilsonTLKernel  apply_from_x(U, x, m );
  

  Ax.apply_everywhere(apply_from_x);
  r0 = b - Ax;

  ScalarFermionField p(x), p0(r0), r(r0), r0star(r), s(x);
  Bulk::WilsonTLKernel apply_from_p(U, p, m );
  Bulk::WilsonTLKernel apply_from_s(U, s, m );

  ScalarFermionField Ap(x), As(x);

  int count = 0;

  nr1 = (r ^ r0star);
  while( count < COUNT_MAX )
    {

      ++count;
      nr = nr1;

      Ap.apply_everywhere(apply_from_p);

      alpha = nr / ( Ap ^ r0star );
      s     = r - Ap * alpha;
      
      As.apply_everywhere(apply_from_s);

      omega = ( As * s ) / ( As * As );
      x     = x + ( (p * alpha) + (s * omega) );
      r     = s - As * omega;
      nr1   = (r * r0star);
      beta  = (nr1 / nr) * (alpha / omega);
      p     = r + (p - Ap * omega) * beta;
      
      std::cout << count << "\t" << abs(x*x) << "\n";
    }

  x = Ap;
}



///////////////////////////////////////////////////
///////////////////////////////////////////////////
///
///  BiCG stab inverter - first attemp. Returns 0 if inversion
///  tolerance is obtained in less than count_max iterations,
///  otherwise stops the inversion and returns 1.
///  1) Is the definition of tolerance correct? 
///  2) Make an Inverter class and this a possible choice of
///  algorithm?
///
///  \author Michele Brambilla <mib.mic@gmail.com>
///  \date Thu Jan 17 11:19:12 2013
template<int DIM, int count_max>
int BiCGstab(  GluonField& U,  ScalarFermionField& src, ScalarFermionField& x,  double& m, const double& tol )
{
  
  x = src;
  Cplx alpha(0.0,0.0), beta, omega(1.0,0.0);
  Cplx rho_1(0.0,0.0), rho_2(1.0,0.0);
  double norm_r, norm_r1;

  //  b.randomize();  // ok, here preconditioning or something else...

  ScalarFermionField b(x), r0(x), Ax(x);

  Bulk::WilsonTLKernel  apply_from_x_blk(U, x, m );
  Lower::WilsonTLKernel apply_from_x_low(U, x, m );
  Upper::WilsonTLKernel apply_from_x_up (U, x, m );
  
  // Ax.apply_everywhere(apply_from_x_blk);
  Ax.apply_on_timeslice(apply_from_x_low, 0);
  //for x_0 != 0 update all directions
  for (int x0 = 1; x0 < T; ++x0)
    for (Direction mu; mu.is_good(); ++mu)
      Ax.apply_on_timeslice(apply_from_x_blk, x0);
  Ax.apply_on_timeslice(apply_from_x_up, T);

  r0 = b - Ax;

//  std::cout << "#\t" << r0*r0 << "\n#----------------\n";
  ScalarFermionField p(r0), r(r0), t(r0), s(x), v(x), rtilde(r0);

  int count = 0;

  while( count < count_max )
    {

      ++count;
      rho_1 = rtilde*r;

      beta = (rho_1/rho_2) * (alpha/omega);

#if DEBUG == 1
      std::cout << "rho_1 = " << rho_1 
		<< ",\tbeta = " << beta
		<< std::endl;
#endif

      p = r + (p - v * omega) * beta;
      Bulk::WilsonTLKernel  apply_from_p_blk(U, p, m );
      Lower::WilsonTLKernel apply_from_p_low(U, p, m );
      Upper::WilsonTLKernel apply_from_p_up (U, p, m );

      //      v.apply_everywhere( apply_from_p_blk );
      v.apply_on_timeslice(apply_from_p_low, 0);
      //for x_0 != 0 update all directions
      for (int x0 = 1; x0 < T; ++x0)
      	for (Direction mu; mu.is_good(); ++mu)
      	  v.apply_on_timeslice(apply_from_p_blk, x0);
      v.apply_on_timeslice(apply_from_p_up, T);

#if DEBUG == 1
      std::cout << "|p-r| = " << (p-r)*(p-r)
		<< std::endl;
#endif

      alpha = rho_1 / (rtilde*v);       

#if DEBUG == 1
      std::cout << "|v| = " << v*v 
		<< ",\talpha = " << alpha
		<< std::endl;
#endif

      s = r - v * alpha;
      Bulk::WilsonTLKernel  apply_from_s_blk(U, s, m );
      Lower::WilsonTLKernel apply_from_s_low(U, s, m );
      Upper::WilsonTLKernel apply_from_s_up (U, s, m );

           // t.apply_everywhere( apply_from_s_blk );
      t.apply_on_timeslice(apply_from_s_low, 0);
      //for x_0 != 0 update all directions
      for (int x0 = 1; x0 < T; ++x0)
      	for (Direction mu; mu.is_good(); ++mu)
      	  t.apply_on_timeslice(apply_from_s_blk, x0);
      t.apply_on_timeslice(apply_from_s_up, T);


      omega = (t*s) / (t*t);

#if DEBUG == 1
      std::cout << "|s| = " << s*s 
		<< "\t,|t| = " << t*t 
		<< ",\tomega = " << omega
		<< std::endl;
#endif

      x += p * alpha + s * omega;
      r  = s - t * omega;

      if (sqrt((r*r).real()) < tol) {
        x += p * alpha;

		// Ax.apply_everywhere(apply_from_x_blk);
	Ax.apply_on_timeslice(apply_from_x_low, 0);
	//for x_0 != 0 update all directions
	for (int x0 = 1; x0 < T; ++x0)
	  for (Direction mu; mu.is_good(); ++mu)
	    Ax.apply_on_timeslice(apply_from_x_blk, x0);
	Ax.apply_on_timeslice(apply_from_x_up, T);

	r0 = b - Ax;

	std::cout << "#\tResidual:" << abs(r0*r0) << "\n#----------------\n";
	return 0;
      }
      
      rho_2 = rho_1;
      
      std::cout << count << "\t|r| = " << abs(r*r) << "\n\n";
    }

  return 1;
}






// template<int DIM, int count_max>
// int CGR(  GluonField& U,  ScalarFermionField& src, ScalarFermionField& x,  double& m, const double& tol )
// {
  
//   x = src;
//   Cplx alpha(0.0,0.0), beta;
//   Cplx norm_rAr, norm_rAr1;

//   //  b.randomize();  // ok, here preconditioning or something else...

//   ScalarFermionField b(x), r0(x), Ax(x);

//   Bulk::WilsonTLKernel  apply_from_x_blk(U, x, m );
//   Lower::WilsonTLKernel apply_from_x_low(U, x, m );
//   Upper::WilsonTLKernel apply_from_x_up (U, x, m );
  
//   // Ax.apply_everywhere(apply_from_x_blk);
//   Ax.apply_on_timeslice(apply_from_x_low, 0);
//   //for x_0 != 0 update all directions
//   for (int x0 = 1; x0 < T; ++x0)
//     for (Direction mu; mu.is_good(); ++mu)
//       Ax.apply_on_timeslice(apply_from_x_blk, x0);
//   Ax.apply_on_timeslice(apply_from_x_up, T);

//   r0 = b - Ax;

// //  std::cout << "#\t" << r0*r0 << "\n#----------------\n";
//   ScalarFermionField p(r0), r(r0), Ar(r0), Ar1(x), Ap(x), rtilde(r0);

//   int count = 0;

//   Bulk::WilsonTLKernel  apply_from_r_blk(U,r,m);
//   Lower::WilsonTLKernel apply_from_r_low(U,r,m);
//   Upper::WilsonTLKernel apply_from_r_up (U,r,m);

//   Ar.apply_on_timeslice(apply_from_r_low, 0);
//   //for x_0 != 0 update all directions
//   for (int x0 = 1; x0 < T; ++x0)
//     for (Direction mu; mu.is_good(); ++mu)
//       Ar.apply_on_timeslice(apply_from_r_blk, x0);
//   Ar.apply_on_timeslice(apply_from_r_up, T);
//   norm_rAr1 = (r*Ar);

//   Bulk::WilsonTLKernel  apply_from_p_blk(U,p,m);
//   Lower::WilsonTLKernel apply_from_p_low(U,p,m);
//   Upper::WilsonTLKernel apply_from_p_up (U,p,m);
//   Ap.apply_on_timeslice(apply_from_p_low, 0);
//   //for x_0 != 0 update all directions
//   for (int x0 = 1; x0 < T; ++x0)
//     for (Direction mu; mu.is_good(); ++mu)
//       Ap.apply_on_timeslice(apply_from_p_blk, x0);
//   Ap.apply_on_timeslice(apply_from_p_up, T);
  
//   while( count < count_max )
//     {

//       ++count;

//       norm_rAr = norm_rAr1;
//       alpha = norm_rAr/(Ap*Ap);

//       x +=  p * alpha;
//       r -= Ap * alpha;

//       Ar.apply_on_timeslice(apply_from_r_low, 0);
//       //for x_0 != 0 update all directions
//       for (int x0 = 1; x0 < T; ++x0)
// 	for (Direction mu; mu.is_good(); ++mu)
// 	  Ar.apply_on_timeslice(apply_from_r_blk, x0);
//       Ar.apply_on_timeslice(apply_from_r_up, T);
      
//       norm_rAr1 = (r*Ar);
//       beta = norm_rAr1/norm_rAr;

//       p = p*beta+ r;
//       Ap = Ap*beta + Ar;


//       if (sqrt((r*r).re) < tol) {
// 	norm_rAr = norm_rAr1;
// 	alpha = norm_rAr/(Ap*Ap);
	
// 	x +=  p * alpha;
	
// 	r0 = b - Ax;

// 	std::cout << "#\tResidual:" << r0*r0 << "\n#----------------\n";
// 	return 0;
//       }
      
//       std::cout << count << "\t|r| = " << r*r << "\n\n";
//     }

//   return 1;
// }







void invertQ( ScalarFermionField& src, GluonField& U, FermionField& dest ) {
  
  std::cout << "#Fermion Operator Inversion" << std::endl;


  double m = .1;

  ScalarFermionField Mxi (src);
  Bulk::WilsonTLKernel  apply_from_src_blk(U, src, m );
  Lower::WilsonTLKernel apply_from_src_low(U, src, m );
  Upper::WilsonTLKernel apply_from_src_up (U, src, m );

  //  Mxi.apply_everywhere(apply_from_src);
  Mxi.apply_on_timeslice(apply_from_src_low, 0);
  //for x_0 != 0 update all directions
  for (int x0 = 1; x0 < T; ++x0)
    for (Direction mu; mu.is_good(); ++mu)
      Mxi.apply_on_timeslice(apply_from_src_blk, x0);
  Mxi.apply_on_timeslice(apply_from_src_up, T);
  
  double tol = 1e-5;

  if( BiCGstab<DIM, COUNT_MAX>( U, Mxi, dest[0], m, tol ) )
    {
      std::cout << "#Error: unable to converge in less than " 
  		<< COUNT_MAX 
  		<< " iterations." 
  		<< std::endl;
    }

}



void test0( ScalarFermionField& src, GluonField& U, ScalarFermionField& dest ) {
  
  double m = .0;

  Bulk::WilsonTLKernel  apply_from_src_blk(U, src, m );
  Lower::WilsonTLKernel apply_from_src_low(U, src, m );
  Upper::WilsonTLKernel apply_from_src_up (U, src, m );

  dest.apply_on_timeslice(apply_from_src_blk, 0);
  //for x_0 != 0 update all directions
  for (int x0 = 1; x0 < T; ++x0)
    for (Direction mu; mu.is_good(); ++mu)
      dest.apply_on_timeslice(apply_from_src_blk, x0);
  dest.apply_on_timeslice(apply_from_src_blk, T);

}

void test1( ScalarFermionField& src, GluonField& U, ScalarFermionField& dest ) {
  
  double m = .1;
  Bulk::WilsonTLKernel  apply_from_src(U, src, m );
  dest.apply_everywhere(apply_from_src);

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
  // for toric boundary: set the time extend to T
  e[0] = T ;

  // fields
  GluonField U(e, 1, 0, nt());
  // for (int t = 0; t < T; ++t){
  //   SetBgfKernel f(t);
  //   U.apply_on_timeslice(f, t);
  // }
  ScalarFermionField Xi  (e, 1, 0, nt());
  FermionField Psi;

  geometry::Geometry<DIM>::raw_pt_t coords={0,0,0,0};
  Point point = Xi.mk_point(coords);

  // source
  int n = 0;
  for ( ScalarFermionField::iterator it = Xi.begin(); it != Xi.end(); ++it, ++n )
    for( Direction mu(0); mu.is_good(); ++mu) 
      (*it)[mu][0] = Cplx(rand()/(double)RAND_MAX,rand()/(double)RAND_MAX);

  // Psi.push_back(ScalarFermionField(e, 1, 0, nt()));
  // test0( Xi, U, Psi[0] );
  // Psi.push_back(ScalarFermionField(e, 1, 0, nt()));
  // test1( Xi, U, Psi[1] );
  // std::cout << Psi[0]*Psi[0] << std::endl;
  // std::cout << Psi[1]*Psi[1] << std::endl;
  // std::cout << (Psi[1] - Psi[0])*(Psi[1] - Psi[0]) << std::endl;

  Psi.push_back(ScalarFermionField(e, 1, 0, nt()));
  invertQ( Xi, U, Psi );
  //std::cout << (Xi - Psi[0])*(Xi - Psi[0]) << std::endl;

  return 0;
}
