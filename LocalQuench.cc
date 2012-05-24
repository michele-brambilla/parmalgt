#include <LocalGluonField.hpp>
#include <Background.h>
#include <Geometry.hpp>
#include <algorithm>
#include <newQCDpt.h>
#include <map>
#include <iostream>
#include <fstream>

const int DIM = 4;
const int ORD = 6;
const int L = 6;
const int T = 6;
const int s = 0;
const int NRUN = 100;
const int MEAS_FREQ = 10;
const double alpha = .05;
const double taug = -.01;
const double stau = 0.1;
const int GF_MODE = 3;

typedef BGptSU3<bgf::AbelianBgf, ORD> ptSU3;
typedef ptt::PtMatrix<ORD> ptsu3;
typedef BGptGluon<bgf::AbelianBgf, ORD, DIM> ptGluon;
typedef pt::Point<DIM> Point;
typedef pt::Direction<DIM> Direction;

typedef fields::LocalField<ptGluon, DIM> GluonField;
typedef GluonField::neighbors_t nt;

struct SetBgf {

  explicit SetBgf(const int& t_in) : t(t_in) { }
  int t;
  void operator()(GluonField& U, const Point& n) const {
    for (Direction mu; mu.is_good(); ++mu)
      U[n][mu].bgf() =  bgf::get_abelian_bgf(t, mu);
  }
};

template <int N>
class LocalGaugeFixing {
public:
  void operator()(GluonField& U, const Point& n) const { 
    do_it(U, n, mode_selektor<N>());
  }
private:
  template <int M> struct mode_selektor { };
  void do_it(GluonField& U, const Point& n, 
             const mode_selektor<1>&) const {
    // exp version
    ptSU3 omega;
    for (Direction mu; mu.is_good(); ++mu)
      omega += U[n][mu] - U[n - mu][mu];
    
    ptSU3 Omega = exp<bgf::AbelianBgf, ORD>( -alpha * omega.reH());
    ptSU3 OmegaDag = exp<bgf::AbelianBgf, ORD>( alpha * omega.reH());
    
    for (Direction mu; mu.is_good(); ++mu){
      U[n][mu] = Omega * U[n][mu];
      U[n - mu][mu] *= OmegaDag;
    }
  }
  void do_it(GluonField& U, const Point& n,
             const mode_selektor<2>&) const {
    ptSU3 omega;
    for (Direction mu; mu.is_good(); ++mu){
      omega += (U[n][mu] * dag(U[n][mu].bgf()) *
            dag( U[n - mu][mu] ) * U[n - mu][mu].bgf());
    }
    ptSU3 Omega = exp<bgf::AbelianBgf, ORD>( alpha * omega.reH());
    ptSU3 OmegaDag = exp<bgf::AbelianBgf, ORD>( -alpha * omega.reH());
    for (Direction mu; mu.is_good(); ++mu){
      U[n][mu] = Omega * U[n][mu];
      U[n - mu][mu] *= OmegaDag;
    }
  }
  void do_it(GluonField& U, const Point& n,
             const mode_selektor<3>&) const {
    ptSU3 omega;
    omega.zero();
    for (Direction mu; mu.is_good(); ++mu){
      ptSU3 Udag = dag(U[n - mu][mu]);
      bgf::AbelianBgf Vdag = U[n][mu].bgf().dag(), V = Udag.bgf().dag();
      omega += U[n][mu]*Vdag*Udag*V;
    }
    ptSU3 Omega = exp<bgf::AbelianBgf, ORD>( -alpha * omega.reH());
    ptSU3 OmegaDag = exp<bgf::AbelianBgf, ORD>( alpha * omega.reH());
    for (Direction mu; mu.is_good(); ++mu){
      U[n][mu] = Omega * U[n][mu];
      U[n - mu][mu] *= OmegaDag;
    }
  }
};


struct LocalGaugeUpdate {
  // for testing, c.f. below
  //static std::vector<MyRand> rands;
  static MyRand Rand;
  Direction mu;
  ptsu3 M; // zero momentum contribution

  explicit LocalGaugeUpdate(const Direction& nu) : mu(nu) { }

  void operator()(GluonField& U, const Point& n) {
    ptSU3 W;
    W.zero();
    for(Direction nu; nu.is_good(); ++nu)
      if(nu != mu)
        W += U[n + mu][nu] *  dag(U[n][nu] * U[n + nu][mu])
          + dag(U[n-nu][mu] * U[n+mu-nu][nu]) * U[n - nu][nu];
    // Close the staple
    W = U[n][mu] * W;
    // DH Feb. 6, 2012
    ptsu3 tmp  = W.reH() * taug; // take to the algebra
    tmp[0] -= stau*SU3rand(Rand); // add noise
    // use this to check if the multithreaded version gives 
    // identical results to the single threaded one
    // tmp[0] -= stau*SU3rand(rands.at(n));
    U[n][mu] = exp<bgf::AbelianBgf, ORD>(tmp)*U[n][mu]; // back to SU3
#pragma omp critical // TODO maybe one can use a reduce or so here
    M += get_q(U[n][mu]); // zero momentum contribution
  }
};

// testing c.f. above
// std::vector<MyRand> LocalGaugeUpdate::rands;

struct MeasureNorm {
  std::vector<double> norm;
  explicit MeasureNorm(const int& order) : norm(order, 0.0) { }
  void operator()(GluonField& U, const Point& n) {
    for (int i = 0; i < norm.size(); ++i)
      norm[i] += U[n].Norm()[i];
  }
};

struct ZeroModeSubtraction {
  Direction mu;
  ptSU3 M; // zero momentum contribution
  ZeroModeSubtraction(const Direction& nu, const ptsu3& N) :
    mu(nu), M(exp<bgf::AbelianBgf, ORD>(-1*reH(N))) { }
  void operator()(GluonField& U, const Point& n) {
    U[n][mu] = M * U[n][mu];
  }
};

MyRand LocalGaugeUpdate::Rand;

struct P_lower {
  ptSU3 val;
  P_lower () : val(bgf::zero<bgf::AbelianBgf>()) { }
  void operator()(GluonField& U, const Point& n){
    Direction t(0);
    for (Direction k(1); k.is_good(); ++k)
      val += U[n][k].bgf() * U[n + k][t] * 
        dag( U[n +t][k] ) * dag( U[n][t] );
  }
};


// This is the plaquette at t = T, arranged such that the derivative
// w.r.t. eta may be inserted at the very end.

struct P_upper {
  ptSU3 val;
  P_upper () : val(bgf::zero<bgf::AbelianBgf>()) { }
  void operator()(GluonField& U, const Point& n){
    Direction t(0);
    ptSU3 tmp;
    for (Direction k(1); k.is_good(); ++k)
      val += dag( U [n + t][k].bgf() ) * dag( U[n][t] ) *
              U[n][k] * U[n + k][t];
  }
};

struct P_spatial {
  ptSU3 val;
  P_spatial () : val(bgf::zero<bgf::AbelianBgf>()) { }
  void operator()(GluonField& U, const Point& n){
    for (Direction k(1); k.is_good(); ++k)
      for (Direction l(k + 1); l.is_good(); ++l)
        val += U[n][k] * U[n + k][l]
          * dag( U[n + l][k] ) * dag( U[n][l] );
  }
};

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

  P_lower Pl; // Plaquettes U_{0,k} at t = 0
  P_upper Pu; // Plaquettes U_{0,k} at t = T - 1
  U.apply_on_slice_with_bnd(Pl, Direction(0), 0);
  U.apply_on_slice_with_bnd(Pu, Direction(0), T-1);
  //meas_on_timeslice(U, 0, Pl);
  //meas_on_timeslice(U, T - 1, Pu);
  P_upper Pm; // Plaquettes U_{0,k} at t = 1, ..., T - 2
  for (int t = 1; t < T - 1; ++t)
    U.apply_on_slice_with_bnd(Pm, Direction(0), t);
    //meas_on_timeslice(U, t, Pm);

  P_spatial Ps; // Spatial plaq. at t = 1, ..., T-1
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


struct Timer {
  double t, tmp;
  static double t_tot;
  Timer () : t(0.0) { };
  void start() { tmp = omp_get_wtime(); };
  void stop() { 
    double elapsed = omp_get_wtime() - tmp; 
    t += elapsed;
    t_tot += elapsed;};
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
  // timing stuff
  typedef std::map<std::string, Timer> tmap;
  tmap timings;
  // these are for testing of the multithreaded version
  //  LocalGaugeUpdate::rands.resize(L*L*L*(T+1));
  //for (int i = 0; i < L*L*L*(T+1); ++i)
  //   LocalGaugeUpdate::rands[i].init(i);
  LocalGaugeUpdate::Rand.init(12312);
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
    SetBgf f(t);
    //U.apply_on_slice_with_bnd(f, Direction(0), t);
    U.apply_on_timeslice(f, t);
  }
  for (int i_ = 0; i_ < NRUN; ++i_){
    if (! (i_ % MEAS_FREQ) ) {
      std::cout << i_;
      MeasureNorm m(ORD + 1);
      U.apply_everywhere(m);
      for (int i = 0 ; i < m.norm.size(); ++i)
        std::cout << " " << m.norm[i];
      std::cout << "\n";
      measure(U);
    }
    ////////////////////////////////////////////////////////
    //
    //  gauge update
    std::vector<LocalGaugeUpdate> gu;
    for (Direction mu; mu.is_good(); ++mu)
      gu.push_back(LocalGaugeUpdate(mu));
    // for x_0 = 0 update the temporal direction only
    //U.apply_on_slice_with_bnd(gu[0], Direction(0), 0);
    timings["Gauge Update"].start();
    U.apply_on_timeslice(gu[0], 0);
    // for x_0 != 0 update all directions
    for (int t = 1; t < T; ++t){
      for (Direction mu; mu.is_good(); ++mu){
        U.apply_on_timeslice(gu[mu], t);
        //U.apply_on_slice_with_bnd(gu[mu], Direction(0), t);
      }
    }
    timings["Gauge Update"].stop();
    ////////////////////////////////////////////////////////
    //
    //  zero mode subtraction
    std::vector<ZeroModeSubtraction> z;
    double vinv = 1./L/L/L/L; // inverse volume
    for (Direction mu; mu.is_good(); ++mu)
      z.push_back(ZeroModeSubtraction(mu, gu[mu].M*vinv));
    //// for x_0 = 0 update the temporal direction only (as above)
    //U.apply_on_slice_with_bnd(z[0], Direction(0), 0);
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
    LocalGaugeFixing<GF_MODE> gf;
    timings["Gauge Fixing"].start();
    for (int t = 1; t < T; ++t)
      U.apply_on_timeslice(gf, t);
      //U.apply_on_slice_with_bnd(gf, Direction(0), t);
    timings["Gauge Fixing"].stop();
  }
  std::cout << "Timings:\n";
  for (tmap::const_iterator i = timings.begin(); i != timings.end();
       ++i){
    pretty_print(i->first, i->second.t);
  }
  pretty_print("TOTAL", Timer::t_tot);
  return 0;
}
