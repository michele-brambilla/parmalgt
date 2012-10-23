#ifndef _KERNELS_H
#define _KERNELS_H
#include <Point.hpp>
#include <LocalField.hpp>
#include <PtTypes.hpp>
#include <newQCDpt.h>
#include <IO.hpp>
#include <uparam.hpp>
#include <Types.h>
#ifdef _OPENMP
#include <omp.h>
#else
namespace kernels {
  int omp_get_max_threads() { return 1; }
  int omp_get_thread_num() { return 0; }
}
#endif
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
///
///  Kernels.
///
///  Kernels can be measurements or updates for lattice fields.
///
///  \author Dirk Hesse <herr.dirk.hesse@gmail.com>
///  \date Thu May 24 17:55:49 2012

namespace kernels {
  
  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  ///
  ///  Kernel for the gauge update.
  ///
  ///  \tparam BGF Background field to use.
  ///  \tparam ORD Perturbative order.
  ///  \tparam DIN Number of space-time dimensions.
  ///
  ///  \author Dirk Hesse <herr.dirk.hesse@gmail.com>
  ///  \date Thu May 24 17:47:43 2012

  template <class BGF, int ORD,int DIM>
  struct GaugeUpdateKernel {
    typedef BGptSU3<BGF, ORD> ptSU3;
    typedef ptt::PtMatrix<ORD> ptsu3;
    typedef BGptGluon<BGF, ORD, DIM> ptGluon;
    typedef pt::Point<DIM> Point;
    typedef pt::Direction<DIM> Direction;
    typedef fields::LocalField<ptGluon, DIM> GluonField;

   
    // for testing, c.f. below
    static std::vector<MyRand> rands;
    //static MyRand Rand;
    Direction mu;
    std::vector<ptsu3> M; // zero momentum contribution
    
    double taug, stau;

    GaugeUpdateKernel(const Direction& nu, const double& t) :
      mu(nu), M(omp_get_max_threads()), taug(t), stau(sqrt(t)) { }
    
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
      ptsu3 tmp  = W.reH() * -taug; // take to the algebra
      //tmp[0] -= stau*SU3rand(Rand); // add noise
      // use this to check if the multithreaded version gives 
      // identical results to the single threaded one
      tmp[0] -= stau*SU3rand(rands.at(n));
      U[n][mu] = exp<bgf::AbelianBgf, ORD>(tmp)*U[n][mu]; // back to SU3
      //#pragma omp critical // TODO maybe one can use a reduce or so here
      //M[omp_get_thread_num()] += get_q(U[n][mu]); 
                                       // zero momentum contribution
    }
    
    //void reduce(){
    //  typename std::vector<ptsu3>::iterator j = M.begin();
    //  if (++j != M.end())
    //    for (; j != M.end(); ++j)
    //      M[0] += *j;
    //}
  };
  template <class C, int N, int M> std::vector<MyRand> 
    kernels::GaugeUpdateKernel<C,N,M>::rands;


  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  ///
  ///  Kernel for gauge fixing.
  ///
  ///  Note I implemented three different methods of which only the
  ///  third one seems to work fine at the moment.
  ///
  ///  \tparam N   Gauge fixing mode. WARNING: Only mode 3 works!
  ///  \tparam BGF Background field to use.
  ///  \tparam ORD Perturbative order.
  ///  \tparam DIN Number of space-time dimensions.
  ///
  ///  \author Dirk Hesse <herr.dirk.hesse@gmail.com>
  ///  \date Thu May 24 17:49:09 2012
  ///

  template <int N, class BGF, int ORD,int DIM>
  class GaugeFixingKernel {
    typedef BGptSU3<BGF, ORD> ptSU3;
    typedef ptt::PtMatrix<ORD> ptsu3;
    typedef BGptGluon<BGF, ORD, DIM> ptGluon;
    typedef pt::Point<DIM> Point;
    typedef pt::Direction<DIM> Direction;
    typedef fields::LocalField<ptGluon, DIM> GluonField;
  public:
    explicit GaugeFixingKernel (const double& a) : alpha (a) { }
    void operator()(GluonField& U, const Point& n) const { 
    do_it(U, n, mode_selektor<N>());
    }
private:
    double alpha;
  template <int M> struct mode_selektor { };
  void do_it(GluonField& U, const Point& n, 
             const mode_selektor<1>&) const {
    // exp version
    ptSU3 omega;
    for (Direction mu; mu.is_good(); ++mu)
      omega += U[n][mu] - U[n - mu][mu];
    
    ptSU3 Omega = exp<bgf::AbelianBgf, ORD>( alpha * omega.reH());
    ptSU3 OmegaDag = exp<bgf::AbelianBgf, ORD>( -alpha * omega.reH());
    
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
  void do_it(GluonField& U, const Point& n, 
             const mode_selektor<4>&) const {
    // exp version
    ptSU3 omega;
    for (Direction mu; mu.is_good(); ++mu)
      omega += U[n][mu]*dag(U[n][mu].bgf()) - 
        dag(U[n-mu][mu].bgf())*U[n - mu][mu];
    
    ptSU3 Omega = exp<bgf::AbelianBgf, ORD>( alpha * omega.reH());
    ptSU3 OmegaDag = exp<bgf::AbelianBgf, ORD>( -alpha * omega.reH());
    
    for (Direction mu; mu.is_good(); ++mu){
      U[n][mu] = Omega * U[n][mu];
      U[n - mu][mu] *= OmegaDag;
    }
  }

};


  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  ///
  ///  Kernel taking care of zero mode subtraction.
  ///
  ///  \tparam BGF Background field to use.
  ///  \tparam ORD Perturbative order.
  ///  \tparam DIN Number of space-time dimensions.
  ///
  ///  \author Dirk Hesse <herr.dirk.hesse@gmail.com>
  ///  \date Thu May 24 17:50:47 2012

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
    //////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////
    ///
    ///  Constructor
    ///
    ///
    ///  \param nu Direction in which the functor will subtract the
    ///  zero modes.
    ///  \param N The sum over the zero modes to be subtracted.
    ///  \author Dirk Hesse <herr.dirk.hesse@gmail.com>
    ///  \date Thu May 24 17:51:17 2012
  ZeroModeSubtractionKernel(const Direction& nu, const ptsu3& N) :
    mu(nu), M(exp<bgf::AbelianBgf, ORD>(-1*reH(N))) { }
  void operator()(GluonField& U, const Point& n) {
    U[n][mu] = M * U[n][mu];
  }
};
  

  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  ///
  ///  Kernel to set the background field to the usual Abelian one.
  ///
  ///  \author Dirk Hesse <herr.dirk.hesse@gmail.com>
  ///  \date Thu May 24 17:52:27 2012
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

  
  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  ///
  ///  Kernel to measure the norm.
  ///
  ///  \author Dirk Hesse <herr.dirk.hesse@gmail.com>
  ///  \date Thu May 24 17:52:56 2012
  template <class BGF, int ORD,int DIM>
  struct MeasureNormKernel {
    typedef BGptSU3<BGF, ORD> ptSU3;
    typedef ptt::PtMatrix<ORD> ptsu3;
    typedef BGptGluon<BGF, ORD, DIM> ptGluon;
    typedef pt::Point<DIM> Point;
    typedef pt::Direction<DIM> Direction;
    typedef fields::LocalField<ptGluon, DIM> GluonField;

    std::vector<typename array_t<double, ORD+1>::Type> norm;
    explicit MeasureNormKernel() : norm(omp_get_max_threads()) { }
    void operator()(GluonField& U, const Point& n) {
      std::vector<double> tmp = U[n].Norm();
      int i = omp_get_thread_num();
      for (int k = 0; k < ORD+1; ++k)
        norm[i][k] += tmp[k];
    }
    void reduce(){
      for (int i = 1, j = omp_get_max_threads(); i < j; ++i)
        for (int k = 0; k < ORD+1; ++k)
          norm[0][k] += norm[i][k];
    }
  };

  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  ///
  ///  Kernel to measure the temporal Plaquette.
  ///
  ///  This measures the temporal plaquette at t = 0, arranged such
  ///   that the derivative w.r.t. eta may be inserted at the very
  ///   end.
  ///
  ///  \author Dirk Hesse <herr.dirk.hesse@gmail.com>
  ///  \date Thu May 24 17:53:07 2012
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
      ptSU3 tmp(bgf::zero<bgf::AbelianBgf>());
      for (Direction k(1); k.is_good(); ++k)
        tmp += U[n][k].bgf() * U[n + k][t] * 
          dag( U[n +t][k] ) * dag( U[n][t] );
#pragma omp critical
      val += tmp;
    }
  
  };
  // measures Udagger * U
  template <class BGF, int ORD,int DIM>
  struct UdagUKernel {

    typedef BGptSU3<BGF, ORD> ptSU3;
    typedef ptt::PtMatrix<ORD> ptsu3;
    typedef BGptGluon<BGF, ORD, DIM> ptGluon;
    typedef pt::Point<DIM> Point;
    typedef pt::Direction<DIM> Direction;
    typedef fields::LocalField<ptGluon, DIM> GluonField;
    ptSU3 val;
    UdagUKernel () : val(bgf::zero<bgf::AbelianBgf>()) { }
    void operator()(GluonField& U, const Point& n){
      ptSU3 tmp(bgf::zero<bgf::AbelianBgf>());
      for (Direction mu; mu.is_good(); ++mu)
        tmp += dag(U[n][mu]) * U[n][mu];
#pragma omp critical
      val += tmp;
    }
  
};
  
  // Measure the temporal paquette
  template <class BGF, int ORD,int DIM>
  struct TemporalPlaqKernel {

    typedef BGptSU3<BGF, ORD> ptSU3;
    typedef ptt::PtMatrix<ORD> ptsu3;
    typedef BGptGluon<BGF, ORD, DIM> ptGluon;
    typedef pt::Point<DIM> Point;
    typedef pt::Direction<DIM> Direction;
    typedef fields::LocalField<ptGluon, DIM> GluonField;
    
    ptSU3 val;
    // \tilde C = - [d_eta C]
    // at the t = 0 side, we have dagger(e^C) and hence an insertion
    // of \tilde C, since C itself is purely imaginary
    BGF Ctilde;

    TemporalPlaqKernel () : val(bgf::zero<bgf::AbelianBgf>()) { }
    
    void operator()(GluonField& U, const Point& n){
      Direction t(0);
      ptSU3 tmp(bgf::zero<bgf::AbelianBgf>());
      for (Direction k(1); k.is_good(); ++k)
        tmp += U[n][t] * U[n + t][k] * 
	  dag( U[n][k] * U[n + k][t] );
#pragma omp critical
      val += tmp;
    }
    
  }; 

  // Kernel to construct the gauge fixing function at t = 0
  template <class BGF, int ORD,int DIM>
  struct GFMeasKernel {

    typedef BGptSU3<BGF, ORD> ptSU3;
    typedef ptt::PtMatrix<ORD> ptsu3;
    typedef BGptGluon<BGF, ORD, DIM> ptGluon;
    typedef pt::Point<DIM> Point;
    typedef pt::Direction<DIM> Direction;
    typedef fields::LocalField<ptGluon, DIM> GluonField;
    
    ptsu3 val;

    void operator()(GluonField& U, const Point& n){
#pragma omp critical
      val += get_q(U[n][Direction(0)]);
    }
    
  };  
  // Kernel to execute the gauge fixing function at t = 0
  // (the one where we used U_0(\vec y, 0) to construct the gf 
  // function)
  template <class BGF, int ORD,int DIM>
  struct GFApplyKernel {

    typedef BGptSU3<BGF, ORD> ptSU3;
    typedef ptt::PtMatrix<ORD> ptsu3;
    typedef BGptGluon<BGF, ORD, DIM> ptGluon;
    typedef pt::Point<DIM> Point;
    typedef pt::Direction<DIM> Direction;
    typedef fields::LocalField<ptGluon, DIM> GluonField;
    
    ptSU3 Omega, OmegaDag;
    
    GFApplyKernel (ptsu3 omega, const double& alpha,
                                 const int& L) { 
      for (int r = 0; r < ORD; ++r)
        for (int i = 0; i < 3; ++i)
          for (int j = 0; j < 3; ++j)
            if (i != j)
              omega[r](i,j) = 0;
      
      omega /= L*L*L;
      Omega = exp<bgf::AbelianBgf, ORD>( -alpha * omega.reH());
      //OmegaDag = exp<bgf::AbelianBgf, ORD>( -alpha * omega.reH());
    }

    void operator()(GluonField& U, const Point& n){
      //for (Direction mu; mu.is_good(); ++mu){
      static Direction t(0);
      U[n][t] = Omega * U[n][t];
        //U[n - mu][mu] *= OmegaDag;
        //}
    }
    
  };
  // Measure the paquette
  template <class BGF, int ORD,int DIM>
  struct PlaqKernel {

    typedef BGptSU3<BGF, ORD> ptSU3;
    typedef ptt::PtMatrix<ORD> ptsu3;
    typedef BGptGluon<BGF, ORD, DIM> ptGluon;
    typedef pt::Point<DIM> Point;
    typedef pt::Direction<DIM> Direction;
    typedef fields::LocalField<ptGluon, DIM> GluonField;
    
    ptSU3 val;
    // \tilde C = - [d_eta C]
    // at the t = 0 side, we have dagger(e^C) and hence an insertion
    // of \tilde C, since C itself is purely imaginary
    BGF Ctilde;

    PlaqKernel () : val(bgf::zero<bgf::AbelianBgf>()) { }
    
    void operator()(GluonField& U, const Point& n){
      ptSU3 tmp(bgf::zero<bgf::AbelianBgf>());
      for (Direction mu; mu.is_good(); ++mu)
	for (Direction nu(mu + 1); nu.is_good(); ++nu)
	  tmp += U[n][mu] * U[n + mu][nu] * 
	    dag( U[n][nu] * U[n + nu][mu] );
#pragma omp critical
      val += tmp;
    }
    
  };

  // Measure Gamma at t = 0
  template <class BGF, int ORD,int DIM>
  struct GammaLowerKernel {

    typedef BGptSU3<BGF, ORD> ptSU3;
    typedef ptt::PtMatrix<ORD> ptsu3;
    typedef BGptGluon<BGF, ORD, DIM> ptGluon;
    typedef pt::Point<DIM> Point;
    typedef pt::Direction<DIM> Direction;
    typedef fields::LocalField<ptGluon, DIM> GluonField;
    
    ptSU3 val;
    // \tilde C = - [d_eta C]
    // at the t = 0 side, we have dagger(e^C) and hence an insertion
    // of \tilde C, since C itself is purely imaginary
    BGF Ctilde;
    explicit GammaLowerKernel (int L) : val(bgf::zero<bgf::AbelianBgf>()) { 
      Cplx ioL(0, 1./L);
      Ctilde[0] = -ioL;
      Ctilde[1] = ioL/2.;
      Ctilde[2] = ioL/2.;
    }
    void operator()(GluonField& U, const Point& n){
      Direction t(0);
      ptSU3 tmp(bgf::zero<bgf::AbelianBgf>());
      for (Direction k(1); k.is_good(); ++k)
        tmp += Ctilde * dag(U[n][k]) * U[n][t] * 
          U[n + t][k] * dag( U[n + k][t] );
#pragma omp critical
      val += tmp;
    }
    
  };
  // Measure Gamma at t = T - a
  template <class BGF, int ORD,int DIM>
  struct GammaUpperKernel {

    typedef BGptSU3<BGF, ORD> ptSU3;
    typedef ptt::PtMatrix<ORD> ptsu3;
    typedef BGptGluon<BGF, ORD, DIM> ptGluon;
    typedef pt::Point<DIM> Point;
    typedef pt::Direction<DIM> Direction;
    typedef fields::LocalField<ptGluon, DIM> GluonField;
    
    ptSU3 val;
    // here, we need [d_eta C'], which is equal to -[d_eta C], hence
    // we can use Ctilde as above
    BGF Ctilde;
    explicit GammaUpperKernel (int L) : val(bgf::zero<bgf::AbelianBgf>()) { 
      Cplx ioL(0, 1./L);
      Ctilde[0] = -ioL;
      Ctilde[1] = ioL/2.;
      Ctilde[2] = ioL/2.;
    }
    void operator()(GluonField& U, const Point& n){
      Direction t(0);
      ptSU3 tmp(bgf::zero<bgf::AbelianBgf>());
      for (Direction k(1); k.is_good(); ++k)
        tmp += Ctilde * U[n + t][k] * dag(U[n + k][t])
          * dag(U[n][k]) * U[n][t];
#pragma omp critical
      val += tmp;
    }
    
};
  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  ///
  ///  Kernel to measure the temporal Plaquette.
  ///
  ///  This measures the plaquette at t = T, arranged such that the
  ///   derivative w.r.t. eta may be inserted at the very end.
  ///
  ///  \author Dirk Hesse <herr.dirk.hesse@gmail.com>
  ///  \date Thu May 24 17:53:07 2012
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
    ptSU3 tmp(bgf::zero<bgf::AbelianBgf>());
    for (Direction k(1); k.is_good(); ++k)
      tmp += dag( U [n + t][k].bgf() ) * dag( U[n][t] ) *
              U[n][k] * U[n + k][t];
#pragma omp cirtical
    val += tmp;
  }
};

  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  ///
  ///  Kernel to measure the spatial Plaquette.
  ///
  ///
  ///  \author Dirk Hesse <herr.dirk.hesse@gmail.com>
  ///  \date Thu May 24 17:53:07 2012
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
      ptSU3 tmp(bgf::zero<bgf::AbelianBgf>());
      for (Direction k(1); k.is_good(); ++k)
        for (Direction l(k + 1); l.is_good(); ++l)
          tmp += U[n][k] * U[n + k][l]
            * dag( U[n + l][k] ) * dag( U[n][l] );
#pragma omp critical
      val += tmp;
    }
};

  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  ///
  ///  Writing a gluon to a file.
  ///
  ///  \author Dirk Hesse <herr.dirk.hesse@gmail.com>
  ///  \date Fri May 25 15:59:06 2012

  template <class BGF, int ORD,int DIM>
  struct FileWriterKernel {

    typedef BGptSU3<BGF, ORD> ptSU3;
    typedef ptt::PtMatrix<ORD> ptsu3;
    typedef BGptGluon<BGF, ORD, DIM> ptGluon;
    typedef pt::Point<DIM> Point;
    typedef pt::Direction<DIM> Direction;
    typedef fields::LocalField<ptGluon, DIM> GluonField;

    explicit FileWriterKernel (uparam::Param& p) : o(p) { }

    void operator()(GluonField& U, const Point& n){
      for (Direction mu(0); mu.is_good(); ++mu)
#pragma omp critical
        U[n][mu].write(o);
    }
    io::CheckedOut o;
  };
  
  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  ///
  ///  Reading a gluon from a file.
  ///
  ///  \author Dirk Hesse <herr.dirk.hesse@gmail.com>
  ///  \date Wed May 30 18:37:03 2012

  template <class BGF, int ORD,int DIM>
  struct FileReaderKernel {

    typedef BGptSU3<BGF, ORD> ptSU3;
    typedef ptt::PtMatrix<ORD> ptsu3;
    typedef BGptGluon<BGF, ORD, DIM> ptGluon;
    typedef pt::Point<DIM> Point;
    typedef pt::Direction<DIM> Direction;
    typedef fields::LocalField<ptGluon, DIM> GluonField;

    typedef GluonField& first_argument_type;
    typedef const Point& second_argument_type;
    typedef void result_type;

    explicit FileReaderKernel (uparam::Param& p) : i(p) { }

    void operator()(GluonField& U, const Point& n){
      for (Direction mu(0); mu.is_good(); ++mu)
        U[n][mu].read(i);
    }
    io::CheckedIn i;
  };

}

#endif
