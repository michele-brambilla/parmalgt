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
#ifdef _OPENMP
#include <omp.h>
#else
#include <time.h>
#endif

// space-time dimensions
const int DIM = 4;
// perturbative order
const int ORD = 6;
// lattice size
int L;// = 6;
int T;// = 6;
// s parameter for staggered
int s;// = 0;
// total number of gauge updates
int NRUN;// = 20;
// frequency of measurements
int MEAS_FREQ;// = 1;
// testing gauge fixing option -- DO NOT TOUCH!
const int GF_MODE = 3;
//
double taug;// = -.01;
double alpha;// = .05;

// some shorthands
typedef bgf::AbelianBgf Bgf_t; // background field
typedef BGptSU3<Bgf_t, ORD> ptSU3; // group variables
typedef ptt::PtMatrix<ORD> ptsu3; // algebra variables
typedef BGptGluon<bgf::AbelianBgf, ORD, DIM> ptGluon; // gluon
typedef pt::Point<DIM> Point;
typedef pt::Direction<DIM> Direction;

// shorthand for gluon field
typedef fields::LocalField<ptGluon, DIM> GluonField;
typedef GluonField::neighbors_t nt;


//
// Make aliases for the Kernels ...
//


// ... for the gauge update/fixing ...
typedef kernels::GaugeUpdateKernel<Bgf_t, ORD, DIM> GaugeUpdateKernel;
// TODO: get rid of this somehow ...
//template <class C, int N, int M> MyRand kernels::GaugeUpdateKernel<C,N, M>::Rand;
template <class C, int N, int M> std::vector<MyRand> kernels::GaugeUpdateKernel<C,N,M>::rands;
typedef kernels::ZeroModeSubtractionKernel<Bgf_t, ORD, DIM> ZeroModeSubtractionKernel;
typedef kernels::GaugeFixingKernel<GF_MODE, Bgf_t, ORD, DIM> GaugeFixingKernel;

// ... to set the background field ...
typedef kernels::SetBgfKernel<Bgf_t, ORD, DIM> SetBgfKernel;


// ... and for the measurements ...
typedef kernels::PlaqLowerKernel<Bgf_t, ORD, DIM> PlaqLowerKernel;
typedef kernels::PlaqUpperKernel<Bgf_t, ORD, DIM> PlaqUpperKernel;
typedef kernels::PlaqSpatialKernel<Bgf_t, ORD, DIM> PlaqSpatialKernel;
typedef kernels::MeasureNormKernel<Bgf_t, ORD, DIM> MeasureNormKernel;

// ... and for the checkpointing.
typedef kernels::FileWriterKernel<Bgf_t, ORD, DIM> FileWriterKernel;
typedef kernels::FileReaderKernel<Bgf_t, ORD, DIM> FileReaderKernel;

// ... some IO
// TODO: Move this as well!
inline void to_bin_file(std::ofstream& of, const Cplx& c){
  of.write(reinterpret_cast<const char*>(&c.re), sizeof(double));
  of.write(reinterpret_cast<const char*>(&c.im), sizeof(double));
}

inline void write_file(const ptSU3& U, const Cplx& tree, const std::string& fname){
  std::ofstream of(fname.c_str(), std::ios_base::app |  std::ios_base::binary);
  to_bin_file(of, tree);
  for (int i = 0; i < ORD; ++i)
    to_bin_file(of, U[i].Tr());
  of.close();
}


// Our measurement
void measure(GluonField &U){
  std::vector<double> nor(ORD*2*3 + 2);
  SU3 dC; // derivative of C w.r.t. eta
  dC(0,0) = Cplx(0, -2./L);
  dC(1,1) = Cplx(0, 1./L);
  dC(2,2) = Cplx(0, 1./L);
  SU3 d_dC; // derivative of C w.r.t. eta and nu
  d_dC(0,0) = Cplx(0, 0);
  d_dC(1,1) = Cplx(0, 2./L);
  d_dC(2,2) = Cplx(0, -2./L);

  // shorthand for V3
  int V3 = L*L*L;

  PlaqLowerKernel Pl; // Plaquettes U_{0,k} at t = 0
  PlaqUpperKernel Pu; // Plaquettes U_{0,k} at t = T - 1
  U.apply_on_slice_with_bnd(Pl, Direction(0), 0);
  U.apply_on_slice_with_bnd(Pu, Direction(0), T-1);
  //meas_on_timeslice(U, 0, Pl);
  //meas_on_timeslice(U, T - 1, Pu);
  PlaqUpperKernel Pm; // Plaquettes U_{0,k} at t = 1, ..., T - 2
  for (int t = 1; t < T - 1; ++t)
    U.apply_on_slice_with_bnd(Pm, Direction(0), t);
    //meas_on_timeslice(U, t, Pm);

  PlaqSpatialKernel Ps; // Spatial plaq. at t = 1, ..., T-1
  for (int t = 1; t < T; ++t)
    U.apply_on_slice_with_bnd(Ps, Direction(0), t);
    //meas_on_timeslice(U, t, Ps);

  // Evaluate Gamma'
  ptSU3 tmp = Pl.val + Pu.val;
  std::for_each(tmp.begin(), tmp.end(), pta::mul(dC));
  Cplx tree = tmp.bgf().ApplyFromLeft(dC).Tr();
  write_file(tmp, tree, "Gp.bindat");

  // Evaluate v
  tmp = Pl.val + Pu.val;
  std::for_each(tmp.begin(), tmp.end(), pta::mul(d_dC));
  tree = tmp.bgf().ApplyFromLeft(d_dC).Tr();
  write_file(tmp, tree, "v.bindat");

  // Evaluate bd_plaq (nomenclature as in MILC code)
  tmp = (Pl.val + Pu.val) / 6 / V3;
  tree = tmp.bgf().Tr();
  write_file(tmp, tree, "bd_plaq.bindat");

  // Evaluate st_plaq
  tmp = (Pl.val + Pu.val + Pm.val) / 3 / V3 / T;
  tree = tmp.bgf().Tr();
  write_file(tmp, tree, "st_plaq.bindat");

  // Eval. ss_plaq
  tmp = Ps.val / 3 / V3 / (T-1);
  tree = tmp.bgf().Tr();
  write_file(tmp, tree, "ss_plaq.bindat");
}


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


// formated cout for the timings
void pretty_print(const std::string& s, const double& d){
  std::cout.width(25); 
  std::cout << s; 
  std::cout.width(0);
  std::cout << ": " << d << "s" << std::endl;
};

double Timer::t_tot = 0.0;

int main(int argc, char *argv[]) {
  // read the parameters
  uparam::Param p;
  p.read("input");
  L = atof(p["L"].c_str());
  s = atoi(p["s"].c_str());
  alpha = atof(p["alpha"].c_str());
  taug = atof(p["taug"].c_str());
  NRUN = atoi(p["NRUN"].c_str());
  MEAS_FREQ = atoi(p["MEAS_FREQ"].c_str());
  T = L;
  // timing stuff
  typedef std::map<std::string, Timer> tmap;
  tmap timings;
  // random stuff
  srand(atoi(p["seed"].c_str()));
  GaugeUpdateKernel::rands.resize(L*L*L*(T+1));
  for (int i = 0; i < L*L*L*(T+1); ++i)
    GaugeUpdateKernel::rands[i].init(rand());
  //GaugeUpdateKernel::Rand.init(12312);
  // generate an array for to store the lattice extents
  geometry::Geometry<DIM>::extents_t e;
  // we want a L = 4 lattice
  std::fill(e.begin(), e.end(), L);
  // for SF boundary: set the time extend to T + 1
  e[0] = T + 1;
  // we will have just one field
  GluonField U(e, 1, 0, nt());
  // initialize background field get method
  bgf::get_abelian_bgf(0, 0, T, L, s);
  // initialzie the background field of U
  for (int t = 0; t <= T; ++t){
    if ( p["read"] == "none" ){
      SetBgfKernel f(t);
      U.apply_on_timeslice(f, t);
    }
    else {
      FileReaderKernel fr(p);
      U.apply_everywhere(fr);
    }
  }
  for (int i_ = 1; i_ <= NRUN; ++i_){
    if (! (i_ % MEAS_FREQ) ) {
      std::cout << i_;
      MeasureNormKernel m(ORD + 1);
      U.apply_everywhere(m);
      for (int i = 0 ; i < m.norm.size(); ++i)
        std::cout << " " << m.norm[i];
      std::cout << "\n";
      timings["measurements"].start();
      measure(U);
      timings["measurements"].stop();
    }
    ////////////////////////////////////////////////////////
    //
    //  gauge update
    std::vector<GaugeUpdateKernel> gu;
    for (Direction mu; mu.is_good(); ++mu)
      gu.push_back(GaugeUpdateKernel(mu, taug));
    timings["Gauge Update"].start();
    // for x_0 = 0 update the temporal direction only
    U.apply_on_timeslice(gu[0], 0);
    // for x_0 != 0 update all directions
    for (int t = 1; t < T; ++t)
      for (Direction mu; mu.is_good(); ++mu)
        U.apply_on_timeslice(gu[mu], t);

    timings["Gauge Update"].stop();
    ////////////////////////////////////////////////////////
    //
    //  zero mode subtraction
    std::vector<ZeroModeSubtractionKernel> z;
    double vinv = 1./L/L/L/L; // inverse volume
    for (Direction mu; mu.is_good(); ++mu){
      gu[mu].reduce();
      z.push_back(ZeroModeSubtractionKernel(mu, gu[mu].M[0]*vinv));
    }
    // for x_0 = 0 update the temporal direction only (as above)
    timings["Zero Mode Subtraction"].start();
    U.apply_on_timeslice(z[0], 0);
    // for x_0 != 0 update all directions (as above)
    for (int t = 1; t < T; ++t)
      for (Direction mu; mu.is_good(); ++mu)
        U.apply_on_timeslice(z[mu], t);
    timings["Zero Mode Subtraction"].stop();

    ////////////////////////////////////////////////////////
    //
    //  gauge fixing
    GaugeFixingKernel gf(alpha);
    timings["Gauge Fixing"].start();
    for (int t = 1; t < T; ++t)
      U.apply_on_timeslice(gf, t);
    timings["Gauge Fixing"].stop();

    //if ( i_ % WR_FREQ == 0 && do_write){
    //  timings["IO"].start();
    //  {
    //    FileWriterKernel fw;
    //    U.apply_everywhere(fw);
    //  }
    //  {
    //    FileReaderKernel fr;
    //    U.apply_everywhere(fr);
    //  }
    //  timings["IO"].stop();
    //}
  } // end main for loop
  if ( p["write"] != "none"){
    FileWriterKernel fw(p);
    U.apply_everywhere(fw);
  }
  // write out timings
  std::cout << "Timings:\n";
  for (tmap::const_iterator i = timings.begin(); i != timings.end();
       ++i){
    pretty_print(i->first, i->second.t);
  }
  pretty_print("TOTAL", Timer::t_tot);
  return 0;
}
