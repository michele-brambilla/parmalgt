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
#include <Clover.hpp>
#include <Methods.hpp>
#ifdef _OPENMP
#include <omp.h>
#else
#include <time.h>
#endif
#include <signal.h>
#include <util.hpp>

#ifdef USE_MPI
#include <mpi.h>
#include <sstream>
#endif

int L;

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
const int ORD = 6;

// some short-hands
typedef bgf::ScalarBgf Bgf_t; // background field
//typedef bgf::AbelianBgf Bgf_t; // background field
typedef BGptSU3<Bgf_t, ORD> ptSU3; // group variables
typedef ptt::PtMatrix<ORD> ptsu3; // algebra variables
typedef BGptGluon<Bgf_t, ORD, DIM> ptGluon; // gluon
typedef pt::Point<DIM> Point; // lattice point
typedef pt::Direction<DIM> Direction; // directions

// shorthand for gluon field
typedef fields::LocalField<ptGluon, DIM> GluonField;
typedef GluonField::neighbors_t nt;

//
// Make aliases for the Kernels ...
//

// ... to set the background field ...
typedef kernels::SetBgfKernel<GluonField> SetBgfKernel;

// ... and for the measurements ...
typedef kernels::MeasureNormKernel<GluonField> MeasureNormKernel;
typedef kernels::PlaqKernel<GluonField> PlaqKernel;

// ... and for the checkpointing.
typedef kernels::FileWriterKernel<GluonField> FileWriterKernel;
typedef kernels::FileReaderKernel<GluonField> FileReaderKernel;

// Our measurement...

// Stuff we want to measure in any case
void measure_common(GluonField &U, const std::string& rep_str){
  // Norm of the Gauge Field
  //MeasureNormKernel m;
  //io::write_file(U.apply_everywhere(m).reduce(),
  //               "SFNorm" + rep_str + ".bindat");
  PlaqKernel p;
  std::list<double> pl, E0, E1;
#ifdef __GXX_EXPERIMENTAL_CXX0X__
  pl.push_back(1);
  for (auto &i : U.apply_everywhere(p).val) pl.push_back(-i.tr().real());
  io::write_file(pl, "SFPlaq" + rep_str + ".bindat");
#endif
  clover::E0m<GluonField> e0m;
  clover::E0s<GluonField> e0s;
  U.apply_on_timeslice(e0m, L/2).reduce();
  U.apply_on_timeslice(e0s, L/2).reduce();
  io::write_file(e0m.result, "E0m.bindat");
  io::write_file(e0s.result, "E0s.bindat");
}

// Stuff that makes sense only for an Abelian background field.
void measure(GluonField &U, const std::string& rep_str,
             const bgf::AbelianBgf&){
  measure_common(U, rep_str);
}

// Stuff that makes sense only for a scalar background field.
void measure(GluonField &U, const std::string& rep_str,
             const bgf::ScalarBgf&){
  measure_common(U, rep_str);
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

double Timer::t_tot = 0.0;

// helper function to feed int, double etc. to the parameters
template <typename T> 
std::string to_string(const T& x){
  std::stringstream sts;
  sts << x;
  return sts.str();
}


int main(int argc, char *argv[]) {
  signal(SIGUSR1, kill_handler);
  signal(SIGUSR2, kill_handler);
  signal(SIGXCPU, kill_handler);
  signal(SIGINT, kill_handler);
  int rank;
#ifdef USE_MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  std::string rank_str = "." + to_string(rank);
#else
  std::string rank_str = "";
#endif
  ////////////////////////////////////////////////////////////////////
  // read the parameters
  uparam::Param p;
  p.read("flow_ipt" + rank_str);
  std::ofstream of(("SFrun.info"+rank_str).c_str(), std::ios::app);
  of << "INPUT PARAMETERS:\n";
  p.print(of);
  of.close();
  L = atoi(p["L"].c_str());  // Spatial lattice size
  int s = atoi(p["s"].c_str());  // s parameter, T = L - s
  double alpha = atof(p["alpha"].c_str());  // gauge fixing parameter
  double taug = atof(p["taug"].c_str()); // integration step size
  double tauf = atof(p["tauf"].c_str());
  int NRUN = atoi(p["NRUN"].c_str()); // Total # of gauge updates
  int NFLOW = atoi(p["NFLOW"].c_str()); // # of wilson flow steps
  int FLOW_FREQ = atoi(p["FLOW_FREQ"].c_str()); // after how many steps shall we flow?
  int MEAS_FREQ = atoi(p["MEAS_FREQ"].c_str()); // freq. of meas.
  int T = L-s; // temporal lattice size
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
  // we want a T x L^3 lattice
  std::fill(e.begin(), e.end(), L);
  // for SF boundary: set the time extend to T + 1
  e[0] = T + 1;
  // initialize background field get method
  bgf::get_abelian_bgf(0, 0, T, L, s);
  // we will have just one field
  GluonField U(e, 1, 0, nt());
  ////////////////////////////////////////////////////////////////////
  //
  // initialzie the background field of U or read config
  if ( p["read"] == "none" )
    for (int t = 0; t <= T; ++t){
      SetBgfKernel f(t);
      U.apply_on_timeslice(f, t);
    }
  else {
    FileReaderKernel fr(p);
    U.apply_everywhere(fr);
  }

  ////////////////////////////////////////////////////////////////////
  //
  // start the simulation
  int i_;
  for (i_ = 1; i_ <= NRUN && !soft_kill; ++i_){
    ////////////////////////////////////////////////////////
    //
    //  gauge update
    ////
    timings["Gauge Update"].start();
    meth::gu::RK2_update(U, taug);
    timings["Gauge Update"].stop();
    ////////////////////////////////////////////////////////
    //
    //  gauge fixing
    timings["Gauge Fixing"].start();
    meth::gf::sf_gauge_fixing(U, alpha);
    timings["Gauge Fixing"].stop();

    ////////////////////////////////////////////////////////
    //
    //  Wilson Flow
    if (! (i_ % FLOW_FREQ) ) {
      // make a copy of the gluon field
      GluonField Up(U);
      // do the flow
      // ... prepare Kernels
      for (int j_ = 1; j_ <= NFLOW && !soft_kill; ++j_){
        timings["Wilson flow"].start();
	meth::wflow::RK2_flow(Up, tauf);
	timings["Wilson flow"].stop();
	if ( ! (j_ % MEAS_FREQ) ){
	  timings["measurements"].start();
	  measure(Up, rank_str, Bgf_t());
	  timings["measurements"].stop();
	}
      }
    }
    
  } // end main for loop
  // write the gauge configuration
  if ( p["write"] != "none"){
    FileWriterKernel fw(p);
    U.apply_everywhere(fw);
  }
  // write out timings
  of.open(("SFrun.info"+rank_str).c_str(), std::ios::app);
  of << "Timings:\n";
  for (tmap::const_iterator i = timings.begin(); i != timings.end();
       ++i){
    util::pretty_print(i->first, i->second.t, "s", of);
  }
  util::pretty_print("TOTAL", Timer::t_tot, "s", of);
  if (soft_kill)
    util::pretty_print("actual # of configs", i_, "", of);
  of.close();

  return 0;
}
