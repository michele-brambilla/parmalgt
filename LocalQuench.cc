#include <LocalGluonField.hpp>
#include <Background.h>
#include <Geometry.hpp>
#include <algorithm>
#include <newQCDpt.h>
#include <iostream>

const int DIM = 4;
const int ORD = 6;
const int L = 4;
const int T = 4;
const int s = 0;
const int NRUN = 1000;
const double alpha = .05;
const double taug = -.01;
const double stau = .1;
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

  static MyRand Rand;
  Direction mu;
  ptsu3 M; // zero momentum contribution
  int count;

  explicit LocalGaugeUpdate(const Direction& nu) : mu(nu), count(0) { }

  void operator()(GluonField& U, const Point& n) {
    ++count;
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
    U[n][mu] = exp<bgf::AbelianBgf, ORD>(tmp)*U[n][mu]; // back to SU3
    M += get_q(U[n][mu]); // zero momentum contribution
  }
};

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

int main(int argc, char *argv[]) {
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
    U.apply_on_slice_with_bnd(f, Direction(0), t);
  }
  for (int i_ = 0; i_ < NRUN; ++i_){
    if (! (i_ % 30) ) {
      MeasureNorm m(ORD + 1);
      U.apply_everywhere(m);
      for (int i = 0 ; i < m.norm.size(); ++i)
        std::cout << m.norm[i] << " ";
      std::cout << "\n";
    }
    ////////////////////////////////////////////////////////
    //
    //  gauge update
    std::vector<LocalGaugeUpdate> gu;
    for (Direction mu; mu.is_good(); ++mu)
      gu.push_back(LocalGaugeUpdate(mu));
    // for x_0 = 0 update the temporal direction only
    U.apply_on_slice_with_bnd(gu[0], Direction(0), 0);
    // for x_0 != 0 update all directions
    for (int t = 1; t < T; ++t){
      for (Direction mu; mu.is_good(); ++mu){
        U.apply_on_slice_with_bnd(gu[mu], Direction(0), t);
      }
    }
    ////////////////////////////////////////////////////////
    //
    //  zero mode subtraction
    std::vector<ZeroModeSubtraction> z;
    double vinv_eff = 1.; // effective inverse volume
    for (Direction mu; mu.is_good(); ++mu) vinv_eff /= gu[mu].count;
    for (Direction mu; mu.is_good(); ++mu)
      z.push_back(ZeroModeSubtraction(mu, gu[mu].M*vinv_eff));
    //// for x_0 = 0 update the temporal direction only (as above)
    U.apply_on_slice_with_bnd(z[0], Direction(0), 0);
    // for x_0 != 0 update all directions (as above)
    for (int t = 1; t < T; ++t)
      for (Direction mu; mu.is_good(); ++mu)
        U.apply_on_slice_with_bnd(z[mu], Direction(0), t);
    ////////////////////////////////////////////////////////
    //
    //  gauge fixing
    LocalGaugeFixing<GF_MODE> gf;
    for (int t = 1; t < T; ++t)
      U.apply_on_slice_with_bnd(gf, Direction(0), t);
    
  }
  return 0;
}
