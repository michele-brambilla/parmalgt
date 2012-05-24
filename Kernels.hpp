#ifndef _KERNELS_H
#define _KERNELS_H
#include <Point.hpp>
#include <LocalField.hpp>
#include <PtTypes.hpp>
#include <newQCDpt.h>

namespace kernels {
  
  // here for now, should move to parameters
  const double taug = -.01;
  const double stau = 0.1;
  const double alpha = .05;

  // singleton Rand

  template <class BGF, int ORD,int DIM>
  struct GaugeUpdateKernel {
    typedef BGptSU3<BGF, ORD> ptSU3;
    typedef ptt::PtMatrix<ORD> ptsu3;
    typedef BGptGluon<BGF, ORD, DIM> ptGluon;
    typedef pt::Point<DIM> Point;
    typedef pt::Direction<DIM> Direction;
    typedef fields::LocalField<ptGluon, DIM> GluonField;

   
    // for testing, c.f. below
    //static std::vector<MyRand> rands;
    static MyRand Rand;
    Direction mu;
    ptsu3 M; // zero momentum contribution
    
    explicit GaugeUpdateKernel(const Direction& nu) : mu(nu) { }
    
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

  template <int N, class BGF, int ORD,int DIM>
  class GaugeFixingKernel {
    typedef BGptSU3<BGF, ORD> ptSU3;
    typedef ptt::PtMatrix<ORD> ptsu3;
    typedef BGptGluon<BGF, ORD, DIM> ptGluon;
    typedef pt::Point<DIM> Point;
    typedef pt::Direction<DIM> Direction;
    typedef fields::LocalField<ptGluon, DIM> GluonField;
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


  template <class BGF, int ORD,int DIM>
  struct ZeroModeSubtractionKernel
  {
    typedef BGptSU3<BGF, ORD> ptSU3;
    typedef ptt::PtMatrix<ORD> ptsu3;
    typedef BGptGluon<BGF, ORD, DIM> ptGluon;
    typedef pt::Point<DIM> Point;
    typedef pt::Direction<DIM> Direction;
    typedef fields::LocalField<ptGluon, DIM> GluonField;
  Direction mu;
  ptSU3 M; // zero momentum contribution
  ZeroModeSubtractionKernel(const Direction& nu, const ptsu3& N) :
    mu(nu), M(exp<bgf::AbelianBgf, ORD>(-1*reH(N))) { }
  void operator()(GluonField& U, const Point& n) {
    U[n][mu] = M * U[n][mu];
  }
};
  

  template <class BGF, int ORD,int DIM>
  struct SetBgfKernel {
    typedef BGptSU3<BGF, ORD> ptSU3;
    typedef ptt::PtMatrix<ORD> ptsu3;
    typedef BGptGluon<BGF, ORD, DIM> ptGluon;
    typedef pt::Point<DIM> Point;
    typedef pt::Direction<DIM> Direction;
    typedef fields::LocalField<ptGluon, DIM> GluonField;

    explicit SetBgfKernel(const int& t_in) : t(t_in) { }
    int t;
    void operator()(GluonField& U, const Point& n) const {
      for (Direction mu; mu.is_good(); ++mu)
        U[n][mu].bgf() =  bgf::get_abelian_bgf(t, mu);
    }
  };


  template <class BGF, int ORD,int DIM>
  struct MeasureNormKernel {
    typedef BGptSU3<BGF, ORD> ptSU3;
    typedef ptt::PtMatrix<ORD> ptsu3;
    typedef BGptGluon<BGF, ORD, DIM> ptGluon;
    typedef pt::Point<DIM> Point;
    typedef pt::Direction<DIM> Direction;
    typedef fields::LocalField<ptGluon, DIM> GluonField;

    std::vector<double> norm;
    explicit MeasureNormKernel(const int& order) : norm(order, 0.0) { }
    void operator()(GluonField& U, const Point& n) {
      for (int i = 0; i < norm.size(); ++i)
        norm[i] += U[n].Norm()[i];
    }
  };

  template <class BGF, int ORD,int DIM>
  struct PlaqLowerKernel {

    typedef BGptSU3<BGF, ORD> ptSU3;
    typedef ptt::PtMatrix<ORD> ptsu3;
    typedef BGptGluon<BGF, ORD, DIM> ptGluon;
    typedef pt::Point<DIM> Point;
    typedef pt::Direction<DIM> Direction;
    typedef fields::LocalField<ptGluon, DIM> GluonField;
    ptSU3 val;
    PlaqLowerKernel () : val(bgf::zero<bgf::AbelianBgf>()) { }
  void operator()(GluonField& U, const Point& n){
    Direction t(0);
    for (Direction k(1); k.is_good(); ++k)
      val += U[n][k].bgf() * U[n + k][t] * 
        dag( U[n +t][k] ) * dag( U[n][t] );
  }
};


// This is the plaquette at t = T, arranged such that the derivative
// w.r.t. eta may be inserted at the very end.

  template <class BGF, int ORD,int DIM>
  struct PlaqUpperKernel {

    typedef BGptSU3<BGF, ORD> ptSU3;
    typedef ptt::PtMatrix<ORD> ptsu3;
    typedef BGptGluon<BGF, ORD, DIM> ptGluon;
    typedef pt::Point<DIM> Point;
    typedef pt::Direction<DIM> Direction;
    typedef fields::LocalField<ptGluon, DIM> GluonField;

  ptSU3 val;
  PlaqUpperKernel () : val(bgf::zero<bgf::AbelianBgf>()) { }
  void operator()(GluonField& U, const Point& n){
    Direction t(0);
    ptSU3 tmp;
    for (Direction k(1); k.is_good(); ++k)
      val += dag( U [n + t][k].bgf() ) * dag( U[n][t] ) *
              U[n][k] * U[n + k][t];
  }
};

  template <class BGF, int ORD,int DIM>
  struct PlaqSpatialKernel {

    typedef BGptSU3<BGF, ORD> ptSU3;
    typedef ptt::PtMatrix<ORD> ptsu3;
    typedef BGptGluon<BGF, ORD, DIM> ptGluon;
    typedef pt::Point<DIM> Point;
    typedef pt::Direction<DIM> Direction;
    typedef fields::LocalField<ptGluon, DIM> GluonField;
  ptSU3 val;
  PlaqSpatialKernel () : val(bgf::zero<bgf::AbelianBgf>()) { }
  void operator()(GluonField& U, const Point& n){
    for (Direction k(1); k.is_good(); ++k)
      for (Direction l(k + 1); l.is_good(); ++l)
        val += U[n][k] * U[n + k][l]
          * dag( U[n + l][k] ) * dag( U[n][l] );
  }
};

}

#endif
