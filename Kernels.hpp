#ifndef _KERNELS_H
#define _KERNELS_H
#include <Point.hpp>
#include <LocalField.hpp>
#include <PtTypes.hpp>
#include <newQCDpt.h>
#include "newMyQCD.h"

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
///  Gamma Matrices.
///
///  Indices and values of all the Dirac gamma and their products
///
///  \author Michele Brambilla <mib.mic@gmail.com>
///  \date Fri Nov 02 12:06:15 2012
///  

namespace gamma{

  int gmuind[15][4]  = {
    {3,2,1,0}, //gm1
    {3,2,1,0},
    {2,3,0,1},
    {2,3,0,1},
    {0,1,2,3}, //gm5
    {3,2,1,0}, //gm51
    {3,2,1,0}, 
    {2,3,0,1},
    {2,3,0,1},
    {0,1,2,3}, // gm12
    {1,0,3,2},
    {1,0,3,2},
    {1,0,3,2},
    {1,0,3,2},
    {0,1,2,3},
  };

  Cplx gmuval[15][4] = {
    {Cplx(0,-1),Cplx(0,-1),Cplx(0, 1),Cplx(0, 1)}, //gm1
    {       -1 ,        1 ,        1 ,       -1 },
    {Cplx(0,-1),Cplx(0, 1),Cplx(0, 1),Cplx(0,-1)},
    {        1 ,        1 ,        1 ,        1 },
    {        1 ,        1 ,       -1 ,       -1 }, // gm5
    {Cplx(0, 1),Cplx(0, 1),Cplx(0, 1),Cplx(0, 1)}, //gm51
    {        1 ,       -1 ,        1 ,       -1 },
    {Cplx(0,1),Cplx(0, -1),Cplx(0, 1),Cplx(0,-1)},
    {       -1 ,       -1 ,        1 ,        1 },
    {Cplx(0,1),Cplx(0, -1),Cplx(0, 1),Cplx(0,-1)}, //gm12
    {       -1 ,        1 ,       -1 ,        1 },
    {Cplx(0,-1),Cplx(0,-1),Cplx(0, 1),Cplx(0, 1)},
    {Cplx(0, 1),Cplx(0, 1),Cplx(0, 1),Cplx(0, 1)},
    {       -1 ,        1 ,        1 ,       -1 },
    {Cplx(0,-1),Cplx(0, 1),Cplx(0, 1),Cplx(0,-1)}
  };


}


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
  
  template <class BGF, int ORD,int DIM>
  struct StapleSqKernel {
    typedef BGptSU3<BGF, ORD> ptSU3;
    typedef ptt::PtMatrix<ORD> ptsu3;
    typedef BGptGluon<BGF, ORD, DIM> ptGluon;
    typedef pt::Point<DIM> Point;
    typedef pt::Direction<DIM> Direction;
    typedef fields::LocalField<ptGluon, DIM> GluonField;
    typedef typename array_t<ptSU3, 1>::Type ptsu3_array_t;
    typedef typename array_t<double, 1>::Type weight_array_t;    

    ptsu3_array_t val;
    Direction mu;
    static weight_array_t weights;

    StapleSqKernel(const Direction& nu) : mu(nu) {  }

    void operator()(GluonField& U, const Point& n) {      
      // std::cout << val[0][0] << "\n";
      val[0].zero();
      // std::cout << val[0][0] << "\n";
      for(Direction nu; nu.is_good(); ++nu)
        if(nu != mu)
          val[0] += U[n + mu][nu] *  dag(U[n][nu] * U[n + nu][mu])
            + dag(U[n-nu][mu] * U[n+mu-nu][nu]) * U[n - nu][nu];
      // Close the staple
      val[0] = U[n][mu] * val[0] ;
    }
    
    ptSU3& reduce() { 
      return val[0]; 
    }

  };
  template <class BGF, int ORD,int DIM>
  typename StapleSqKernel<BGF,ORD,DIM>::weight_array_t 
  StapleSqKernel<BGF,ORD,DIM>::weights;

  
  ////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////
  ///
  /// Kernel to measure the staples for e.g. LW or Iwasaki aciton.
  ///
  /// Be careful when you use this to perform a gauge update. As was
  /// pointed out by Aoki et al. in hep-lat/9808007, generally the
  /// e.o.m do not hold if one naively uses the rectangular
  /// plaquettes. For a perturbative calculation, the weights of the
  /// rectangular loops at the boundary should be adjusted according
  /// to choice 'B' in the paper, c.f. eqns. (2.17) and (2.18).
  ///
  /// This is implemented in the pre- and post-processing kernels
  /// LWProcessA and LWProcessB.
  ///
  /// Commenty by D.H., Oct. 26, 2012

  template <class BGF, int ORD,int DIM>
  struct StapleReKernel {
    typedef BGptSU3<BGF, ORD> ptSU3;
    typedef ptt::PtMatrix<ORD> ptsu3;
    typedef BGptGluon<BGF, ORD, DIM> ptGluon;
    typedef pt::Point<DIM> Point;
    typedef pt::Direction<DIM> Direction;
    typedef fields::LocalField<ptGluon, DIM> GluonField;

    typedef typename array_t<ptSU3, 2>::Type ptsu3_array_t;    
    typedef typename array_t<double, 2>::Type weight_array_t;    

    ptsu3_array_t val;
    Direction mu;
    static weight_array_t weights;
    
    StapleReKernel(const Direction& nu) : mu(nu) {  }

    ptSU3 two_by_one(GluonField& U, const Point& n, const Direction& nu){
      return U[n + mu][mu] * U[n + mu + mu][nu] * 
             dag(U[n][nu] * U[n+nu][mu] * U[n+nu+mu][mu])
           + U[n+mu][mu] * dag( U[n-nu][mu] * 
             U[n+mu-nu][mu] * U[n+mu+mu-nu][nu] ) * U[n-nu][nu]
           + U[n+mu][nu] * dag( U[n-mu][nu] * 
             U[n-mu+nu][mu] * U[n+nu][mu] ) * U[n-mu][mu]
           + dag( U[n-nu-mu][mu] * U[n-nu][mu] * U[n-nu+mu][nu] ) * 
             U[n-nu-mu][nu] * U[n-mu][mu];
    }
    ptSU3 one_by_two(GluonField& U, const Point& n, const Direction& nu){
      return U[n + mu][nu] * U[n+mu+nu][nu] * 
             dag(U[n][nu] * U[n+nu][nu] * U[n+nu+nu][mu]) 
           + dag( U[n-nu-nu][mu] * U[n-nu-nu+mu][nu] * U[n-nu+mu][nu] ) * 
             U[n-nu-nu][nu] * U[n-nu][nu];
    }
    ptSU3 one_by_one(GluonField& U, const Point& n, const Direction& nu){
       return U[n + mu][nu] *  dag(U[n][nu] * U[n + nu][mu])
	    + dag(U[n-nu][mu] * U[n+mu-nu][nu]) * U[n - nu][nu];
    }
    
    void operator()(GluonField& U, const Point& n) {      

      // is it better to use a foreach ?
      for( int i = 0; i < val.size(); ++i) {
	val[i].bgf() *= 0;
	val[i].zero();
      }

      for(Direction nu; nu.is_good(); ++nu)
        if(nu != mu) {
	  // 1x1 contribution
	  val[0] += one_by_one(U, n, nu);
	  // 2x1 contribution
	  val[1] += two_by_one(U, n, nu) + one_by_two(U, n, nu);
	}
      
      // Close the staple
      for( int i = 0; i < val.size(); ++i) {
	val[i] = U[n][mu] * val[i] ;
      }

    }
    ptSU3& reduce() { 
      val[0] *= weights[0];
      for( int i = 1; i < val.size(); ++i) val[0] += weights[i]*val[i] ;
      return val[0]; 
    }

  };

  template <class BGF, int ORD,int DIM>
  typename StapleReKernel<BGF,ORD,DIM>::weight_array_t StapleReKernel<BGF,ORD,DIM>::weights;

  template <class BGF, int ORD,int DIM>
  struct TrivialPreProcess {
    typedef BGptSU3<BGF, ORD> ptSU3;
    typedef ptt::PtMatrix<ORD> ptsu3;
    typedef BGptGluon<BGF, ORD, DIM> ptGluon;
    typedef pt::Point<DIM> Point;
    typedef pt::Direction<DIM> Direction;
    typedef fields::LocalField<ptGluon, DIM> GluonField;
    static void pre_process (GluonField& , const Point&, const Direction&) { }
    static void post_process (GluonField& , const Point&, const Direction& ) { }
  };

  // To be used at x with x_0 = a.
  // This modifies the background field, such that O(a) improvement
  // and the equations of motion hold at tree level.
  template <class BGF, int ORD,int DIM>
  struct LWProcessA {
    typedef BGptSU3<BGF, ORD> ptSU3;
    typedef ptt::PtMatrix<ORD> ptsu3;
    typedef BGptGluon<BGF, ORD, DIM> ptGluon;
    typedef pt::Point<DIM> Point;
    typedef pt::Direction<DIM> Direction;
    typedef fields::LocalField<ptGluon, DIM> GluonField;

    static void pre_process (GluonField& U, const Point& n, const Direction& k) { 
      static Direction t(0);
      U[n - t + k + k][t] *= 1.5;
      U[n - t - k][t] *= 1.5;
    }
    static void post_process (GluonField& U, const Point& n, const Direction& k) { 
      static Direction t(0);
      U[n - t + k + k][t] /= 1.5;
      U[n - t - k][t] /= 1.5;
    }
  };

  // To be used at x with x_0 = T - a.
  // This modifies the background field, such that O(a) improvement
  // and the equations of motion hold at tree level.
  template
  <class BGF, int ORD,int DIM>
  struct LWProcessB {
    typedef BGptSU3<BGF, ORD> ptSU3;
    typedef ptt::PtMatrix<ORD> ptsu3;
    typedef BGptGluon<BGF, ORD, DIM> ptGluon;
    typedef pt::Point<DIM> Point;
    typedef pt::Direction<DIM> Direction;
    typedef fields::LocalField<ptGluon, DIM> GluonField;
    
    static void pre_process (GluonField& U, const Point& n, const Direction& k) { 
      static Direction t(0);
      U[n - k][t] *= 1.5; // ...
      U[n + k + k][t] *= 1.5;
    }
    static void post_process (GluonField& U, const Point& n, const Direction& k) { 
      static Direction t(0);
      U[n - k][t] /= 1.5; // ...
      U[n + k + k][t] /= 1.5;
    }
  };


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

  template <class BGF, int ORD,int DIM, 
	    class StapleK_t, class Process >
  struct GaugeUpdateKernel {

    typedef BGptSU3<BGF, ORD> ptSU3;
    typedef ptt::PtMatrix<ORD> ptsu3;
    typedef BGptGluon<BGF, ORD, DIM> ptGluon;
    typedef pt::Point<DIM> Point;
    typedef pt::Direction<DIM> Direction;
    typedef fields::LocalField<ptGluon, DIM> GluonField;

    typedef std::vector<Cplx>::iterator cpx_vec_it;
    typedef std::vector<std::vector<Cplx> >::iterator outer_cvec_it;

    // for testing, c.f. below
    static std::vector<MyRand> rands;
    //static MyRand Rand;
    Direction mu;

#ifndef SF
    std::vector<ptsu3> M; // zero momentum contribution
#endif
    
    double taug, stau;

    // on-the-fly plaquette measure
    std::vector<std::vector<Cplx> > plaq, pp;

    GaugeUpdateKernel(const Direction& nu, const double& t) :
      mu(nu), 
#ifndef SF
      M(omp_get_max_threads()), 
#endif
      taug(t), stau(sqrt(t)), 
      plaq(omp_get_max_threads(), std::vector<Cplx>(ORD+1)), 
      pp(omp_get_max_threads(), std::vector<Cplx>(ORD+1)) { }

    void operator()(GluonField& U, const Point& n) {
      ptSU3 W;

      // We wants this static, but it fails ... field grows bigger and bigger ...

      // Make a Kernel to calculate and store the plaquette(s)
      StapleK_t st(mu); // maye make a vector of this a class member
      Process::pre_process(U,n,mu);
      st(U,n);
      Process::post_process(U,n,mu);
      pp[omp_get_thread_num()] = (st.val[0].trace()); // Save the 1x1 plaquette

      for(cpx_vec_it k = plaq[omp_get_thread_num()].begin(), 
	    j = pp[omp_get_thread_num()].begin(),
	    e = plaq[omp_get_thread_num()].end(); 
	  k != e; ++k, ++j) *k += *j;

      ptsu3 tmp  = st.reduce().reH() * -taug;

      // DH Feb. 6, 2012
      // ptsu3 tmp  = W.reH() * -taug; // take to the algebra
      //tmp[0] -= stau*SU3rand(Rand); // add noise
      // use this to check if the multithreaded version gives 
      // identical results to the single threaded one
      tmp[0] -= stau*SU3rand(rands.at(n));
      U[n][mu] = exp<BGF, ORD>(tmp)*U[n][mu]; // back to SU3
      //#pragma omp critical // TODO maybe one can use a reduce or so here
#ifndef SF
      M[omp_get_thread_num()] += get_q(U[n][mu]); // zero momentum contribution
#endif
    }
    
    void reduce(){
      outer_cvec_it i = plaq.begin();
      outer_cvec_it pe = plaq.end();
      if (++i != plaq.end())
	for (; i != pe; ++i)
	  plaq[0] += *i;
#ifndef SF
      typename std::vector<ptsu3>::iterator j = M.begin();
      if (++j != M.end())
        for (; j != M.end(); ++j)
          M[0] += *j;
#endif
    }
    
  };
  template <class C, int N, int M, class P, class Q> std::vector<MyRand> 
  kernels::GaugeUpdateKernel<C,N,M,P,Q>::rands;


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
    
    ptSU3 Omega = exp<BGF, ORD>( alpha * omega.reH());
    ptSU3 OmegaDag = exp<BGF, ORD>( -alpha * omega.reH());
    
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
    ptSU3 Omega = exp<BGF, ORD>( alpha * omega.reH());
    ptSU3 OmegaDag = exp<BGF, ORD>( -alpha * omega.reH());
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
      BGF Vdag = U[n][mu].bgf().dag(), V = Udag.bgf().dag();
      omega += U[n][mu]*Vdag*Udag*V;
    }
    ptSU3 Omega = exp<BGF, ORD>( -alpha * omega.reH());
    ptSU3 OmegaDag = exp<BGF, ORD>( alpha * omega.reH());
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
    
    ptSU3 Omega = exp<BGF, ORD>( alpha * omega.reH());
    ptSU3 OmegaDag = exp<BGF, ORD>( -alpha * omega.reH());
    
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
    mu(nu), M(exp<BGF, ORD>(-1*reH(N))) { }
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


  // MiB
  template <int ORD,int DIM>
  struct SetBgfKernel<bgf::ScalarBgf, ORD, DIM> {

    typedef BGptGluon<bgf::ScalarBgf, ORD, DIM> ptGluon;
    typedef pt::Point<DIM> Point;
    typedef pt::Direction<DIM> Direction;
    typedef fields::LocalField<ptGluon, DIM> GluonField;

    explicit SetBgfKernel(const int& t_in) { }

    void operator()(GluonField & U, const Point& n) const { }

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
    PlaqLowerKernel () : val(bgf::zero<BGF>()) { }
    void operator()(GluonField& U, const Point& n){
      Direction t(0);
      ptSU3 tmp(bgf::zero<BGF>());
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
    UdagUKernel () : val(bgf::zero<BGF>()) { }
    void operator()(GluonField& U, const Point& n){
      ptSU3 tmp(bgf::zero<BGF>());
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

    TemporalPlaqKernel () : val(bgf::zero<BGF>()) { }
    
    void operator()(GluonField& U, const Point& n){
      Direction t(0);
      ptSU3 tmp(bgf::zero<BGF>());
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
      Omega = exp<BGF, ORD>( -alpha * omega.reH());
      //OmegaDag = exp<BGF, ORD>( -alpha * omega.reH());
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

    PlaqKernel () : val(bgf::zero<BGF>()) { }
    
    void operator()(GluonField& U, const Point& n){
      ptSU3 tmp(bgf::zero<BGF>());
      for (Direction mu; mu.is_good(); ++mu)
	for (Direction nu(mu + 1); nu.is_good(); ++nu)
	  tmp += U[n][mu] * U[n + mu][nu] * 
	    dag( U[n][nu] * U[n + nu][mu] );
#pragma omp critical
      val += tmp;
    }
    
  };

#ifdef IMP_ACT
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
    explicit GammaLowerKernel (int L) : val(bgf::zero<BGF>()) { 
      Cplx ioL(0, 1./L);
      Ctilde[0] = -ioL;
      Ctilde[1] = ioL/2.;
      Ctilde[2] = ioL/2.;
    }
    void operator()(GluonField& U, const Point& n){
      Direction t(0);
      ptSU3 tmp(bgf::zero<BGF>());
      for (Direction k(1); k.is_good(); ++k){
        // the 1x1 contribution
        tmp += (Ctilde * dag(U[n][k]) * U[n][t] * 
                U[n + t][k] * dag( U[n + k][t] )) * (5./3);
        // the 2x1 contribution w\ two links @ the boundary
        // we actually have two contributions to the derivative form
        // each boundary link. However, we do not care because the
        // insertion of \tilde C commutes with the links. Hence, we
        // just multiply the weight with two, 
        //  w = 2 * (3/2) * c_1 = 3*(-1/12) = -1/4
        tmp += (Ctilde * dag(U[n][k]) * dag(U[n-k][k]) 
                * U[n-k][t] * U[n + t - k][k]
                * U[n + t][k] * dag( U[n + k][t] )) * (-1./4);
        // the 1x2 contribution, with usual weight c_1
        tmp += (Ctilde * dag(U[n][k]) * U[n][t] 
                * U[n+t][t] * U[n+t+t][k] 
                * dag(U[n + k + t][t]) * dag(U[n + k][t])) * (-1./12);
      }
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
    explicit GammaUpperKernel (int L) : val(bgf::zero<BGF>()) { 
      Cplx ioL(0, 1./L);
      Ctilde[0] = -ioL;
      Ctilde[1] = ioL/2.;
      Ctilde[2] = ioL/2.;
    }
    void operator()(GluonField& U, const Point& n){
      Direction t(0);
      ptSU3 tmp(bgf::zero<BGF>());
      for (Direction k(1); k.is_good(); ++k){
        // the 1x1 contribution
        tmp += (Ctilde * U[n + t][k] * dag(U[n + k][t])
                * dag(U[n][k]) * U[n][t]) * (5./3);
        // the 2x1 contribution w\ two links @ the boundary
        // we actually have two contributions to the derivative form
        // each boundary link. However, we do not care because the
        // insertion of \tilde C commutes with the links. Hence, we
        // just multiply the weight with two, 
        //  w = 2 * (3/2) * c_1 = 3*(-1/12) = -1/4
        tmp += (Ctilde * U[n + t][k] * U[n + t + k][k]
                * dag(U[n + k + k][t]) * dag(U[n + k][k])
                * dag(U[n][k]) * U[n][t]) * (1./4);
        // the 1x2 contribution, with usual weight c_1
        tmp += (Ctilde * U[n + t][k] * dag(U[n + k][t])
                * dag(U[n + k - t][t]) * dag(U[n-t][k]) 
                * U[n - t][t] * U[n][t]) * (-1./12);
      }
#pragma omp critical
      val += tmp;
    }
  };


#else
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
    explicit GammaLowerKernel (int L) : val(bgf::zero<BGF>()) { 
      Cplx ioL(0, 1./L);
      Ctilde[0] = -ioL;
      Ctilde[1] = ioL/2.;
      Ctilde[2] = ioL/2.;
    }
    void operator()(GluonField& U, const Point& n){
      Direction t(0);
      ptSU3 tmp(bgf::zero<BGF>());
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
    explicit GammaUpperKernel (int L) : val(bgf::zero<BGF>()) { 
      Cplx ioL(0, 1./L);
      Ctilde[0] = -ioL;
      Ctilde[1] = ioL/2.;
      Ctilde[2] = ioL/2.;
    }
    void operator()(GluonField& U, const Point& n){
      Direction t(0);
      ptSU3 tmp(bgf::zero<BGF>());
      for (Direction k(1); k.is_good(); ++k)
        tmp += Ctilde * U[n + t][k] * dag(U[n + k][t])
          * dag(U[n][k]) * U[n][t];
#pragma omp critical
      val += tmp;
    }
  };
#endif


  // MB: these makes no sense in case of ScalarBgf
  template <int ORD,int DIM>
  struct GammaLowerKernel<bgf::ScalarBgf, ORD, DIM> {
    
    typedef BGptSU3<bgf::ScalarBgf, ORD> ptSU3;
    typedef BGptGluon<bgf::ScalarBgf, ORD, DIM> ptGluon;
    typedef pt::Point<DIM> Point;
    typedef fields::LocalField<ptGluon, DIM> GluonField;
    
    ptSU3 val;

    explicit GammaLowerKernel (int L) { }
    void operator()(GluonField& U, const Point& n) { }
  };

  template <int ORD,int DIM>
  struct GammaUpperKernel<bgf::ScalarBgf, ORD, DIM> {
    
    typedef BGptSU3<bgf::ScalarBgf, ORD> ptSU3;
    typedef BGptGluon<bgf::ScalarBgf, ORD, DIM> ptGluon;
    typedef pt::Point<DIM> Point;
    typedef fields::LocalField<ptGluon, DIM> GluonField;
    
    ptSU3 val;

    explicit GammaUpperKernel (int L) { }
    void operator()(GluonField& U, const Point& n) { }
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
  PlaqUpperKernel () : val(bgf::zero<BGF>()) { }
  void operator()(GluonField& U, const Point& n){
    Direction t(0);
    ptSU3 tmp(bgf::zero<BGF>());
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
    PlaqSpatialKernel () : val(bgf::zero<BGF>()) { }
    void operator()(GluonField& U, const Point& n){
      ptSU3 tmp(bgf::zero<BGF>());
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



  /*
  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  ///
  ///  Perturbative application of the Wilson dirac operator
  ///
  ///  \author Michele Brambilla <mib.mic@gmail.com>
  ///  \date Fri Nov 02 16:24:21 2012


  template <class BGF, int ORD,int DIM>
  struct WilsonPTKernel {

    // geometry
    typedef typename pt::Point<DIM> Point;
    typedef pt::Direction<DIM> Direction;

    // gauge
    typedef BGptSU3<BGF, ORD> ptSU3;
    typedef ptt::PtMatrix<ORD> ptsu3;
    typedef BGptGluon<BGF, ORD, DIM> ptGluon;
    typedef fields::LocalField<ptGluon, DIM> GluonField;

    // fermion
    typedef SpinColor<4> Fermion;
    typedef fields::LocalField<Fermion, DIM> ScalarFermionField;
    typedef std::vector<ScalarFermionField> FermionField;
    typedef typename array_t<double, DIM>::Type double_array_t;    
    typedef typename array_t<double, ORD>::Type ptarray_t;    
    
    WilsonPTKernel(GluonField& G, FermionField& X, ptarray_t& m ) : U(G), src(X), mass(m), ord_(0) { };
    

    void operator() ( ScalarFermionField& dest, Point& n) {

      for( int kord = 0; kord < ord_; ++kord ) {
	// DEBUG
	// if( n == 0 )
	//   std::cout << "F[" << ord_      << "] += "
	// 	    << "U[" << ord_-kord << "]*"
	// 	    << "f[" << kord      << "]"
	// 	    << "\n";

	dest[n] += (src[kord][n] * mass[ord_-kord]);
	Fermion Xi1, Xi2;
	
	// Time direction
	Point dn = n-Direction(0);
	Point up = n+Direction(0);
	// (1 +(-) \gamma_\mu)\psi
	for( Direction nu(0); nu.is_good(); ++nu )
	  {
	    Xi1[nu] = ( src[kord][dn][nu] + gamma::gmuval[0][nu] * src[kord][dn][gamma::gmuind[0][nu]] );
	    Xi2[nu] = ( src[kord][up][nu] - gamma::gmuval[0][nu] * src[kord][up][gamma::gmuind[0][nu]] );
	  }
	// sign takes care of (anti)periodic (or others) boundary conditions
	dest[n] -= ( dag(U[dn][Direction(0)][ord_-kord-1]) * Xi1 * sign[0] + 
		         U[n ][Direction(0)][ord_-kord-1]  * Xi2 * sign[1]) * .5;
	
	// Spatial directions
	for( Direction mu(1); mu.is_good(); ++mu )
	  {
	    
	    Point dn = n-Direction(mu);
	    Point up = n+Direction(mu);
	    // (1 +(-) \gamma_\mu)\psi
	    for( Direction nu(0); nu.is_good(); ++nu )
	      {
		Xi1[nu] = ( src[kord][dn][nu] + gamma::gmuval[mu][nu] * src[kord][dn][gamma::gmuind[mu][nu]] );
		Xi2[nu] = ( src[kord][up][nu] - gamma::gmuval[mu][nu] * src[kord][up][gamma::gmuind[mu][nu]] );
	      }
	    
	    dest[n] -= (  dag(U[dn][Direction(0)][ord_-kord-1]) * Xi1 + 
			      U[n ][Direction(0)][ord_-kord-1]  * Xi2 ) * .5;
	  } // mu

      } // kord

    }



  
    void bulk()  { sign[0] = Cplx( 1,0); sign[1] = Cplx( 1,0); }
    void lower() { sign[0] = Cplx(-1,0); sign[1] = Cplx( 1,0); }
    void upper() { sign[0] = Cplx( 1,0); sign[1] = Cplx(-1,0); }


    inline int& ord() { return ord_; }
    void operator++() { ++ord_; }

  private:
    int ord_;
    ptarray_t mass;
    array_t<Cplx, 2>::Type sign;
    GluonField& U;
    FermionField& src;

  };




  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  ///
  ///  TreeLevel ( U = bgf() ) application of the Wilson dirac operator
  ///
  ///  \author Michele Brambilla <mib.mic@gmail.com>
  ///  \date Fri Nov 02 16:23:58 2012


  template <class BGF, int ORD,int DIM>
  struct WilsonTreeLevelKernel {

    // geometry
    typedef typename pt::Point<DIM> Point;
    typedef pt::Direction<DIM> Direction;

    // gauge
    typedef BGptSU3<BGF, ORD> ptSU3;
    typedef ptt::PtMatrix<ORD> ptsu3;
    typedef BGptGluon<BGF, ORD, DIM> ptGluon;
    typedef fields::LocalField<ptGluon, DIM> GluonField;

    // fermion
    typedef SpinColor<4> Fermion;
    typedef fields::LocalField<Fermion, DIM> ScalarFermionField;
    typedef std::vector<ScalarFermionField> FermionField;
    typedef typename array_t<double, DIM>::Type double_array_t;    
    typedef typename array_t<double, ORD>::Type ptarray_t;    
    
    WilsonTreeLevelKernel(GluonField& G, ScalarFermionField& X, double& m ) : U(G), src(X), mass(m) { };
    
    void operator() ( ScalarFermionField& dest, Point& n) {

      dest[n] = (src[n] * mass);
      Fermion Xi1, Xi2;
      
      // Time direction
      Point dn = n-Direction(0);
      Point up = n+Direction(0);
      // (1 +(-) \gamma_\mu)\psi
      for( Direction nu(0); nu.is_good(); ++nu )
	{
	  Xi1[nu] = ( src[dn][nu] + gamma::gmuval[0][nu] * src[dn][gamma::gmuind[0][nu]] );
	  Xi2[nu] = ( src[up][nu] - gamma::gmuval[0][nu] * src[up][gamma::gmuind[0][nu]] );
	}
      // sign takes care of (anti)periodic (or others) boundary conditions
      dest[n] -= ( dag(U[dn][Direction(0)].bgf()) * Xi1 * sign[0] + 
		       U[n ][Direction(0)].bgf()  * Xi2 * sign[1]) * .5;
      
      // Spatial directions
      for( Direction mu(1); mu.is_good(); ++mu )
	{
	  
	  Point dn = n-Direction(mu);
	  Point up = n+Direction(mu);
	  // (1 +(-) \gamma_\mu)\psi
	  for( Direction nu(0); nu.is_good(); ++nu )
	    {
	      Xi1[nu] =  ( src[dn][nu] + gamma::gmuval[mu][nu] * src[dn][gamma::gmuind[mu][nu]] );
	      Xi2[nu] =  ( src[up][nu] - gamma::gmuval[mu][nu] * src[up][gamma::gmuind[mu][nu]] );
	    }
	  
	  dest[n] -= ( dag(U[dn][mu].bgf()) * Xi1 + 
		           U[n ][mu].bgf()  * Xi2 ) * .5;
	} // mu
      
    }
      
    void bulk()  { sign[0] = Cplx( 1,0); sign[1] = Cplx( 1,0); }
    void lower() { sign[0] = Cplx(-1,0); sign[1] = Cplx( 1,0); }
    void upper() { sign[0] = Cplx( 1,0); sign[1] = Cplx(-1,0); }


    inline int& ord() { return kord; }
    
  private:
    int kord;
    double mass;
    array_t<Cplx, 2>::Type sign;
    GluonField& U;
    ScalarFermionField& src;

  };




  
  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  ///
  ///  Perturbative application of the staggered dirac operator
  ///  
  ///
  ///  \author Michele Brambilla <mib.mic@gmail.com>
  ///  \date Fri Nov 02 16:23:46 2012


  template <class BGF, int ORD,int DIM>
  struct StaggeredPTKernel {

    // geometry
    typedef pt::Point<DIM> Point;
    typedef pt::Direction<DIM> Direction;

    // gauge
    typedef BGptSU3<BGF, ORD> ptSU3;
    typedef ptt::PtMatrix<ORD> ptsu3;
    typedef BGptGluon<BGF, ORD, DIM> ptGluon;
    typedef fields::LocalField<ptGluon, DIM> GluonField;

    // fermion
    typedef SpinColor<1> Fermion;
    typedef fields::LocalField<Fermion, DIM> FermionField;
    typedef typename array_t<double, DIM>::Type double_array_t;    
    
    StaggeredPTKernel( GluonField& G, FermionField& X, const double& m ) : U(G), src(X), mass(m) { Uorder = 0; std::cout << "Fix eta definition!!!!!\n";};
    
    void eta(const Point& n, double_array_t& phase) {

    }
    
    void next_Uorder() { Uorder++; }
    void reset_Uorder() { Uorder = 0; }

    void operator()( FermionField& dest, Point& n) {
      eta(n, phase);
      for( Direction nu(0); nu.is_good(); ++nu )
	dest[n][0] += .5*phase[(int)nu] * ( U[n][nu][Uorder] * src[n+nu][0] - U[n-nu][nu][Uorder] * src[n-nu][0] );
      //      mass
      dest[n][0] += mass * src[n][0];

    }

  private:
    int Uorder;
    double mass;
    double_array_t phase;
    GluonField& U;
    FermionField& src;

  };


  /*
  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  ///
  ///  Tree-level application of the staggered dirac operator
  ///
  ///  \author Michele Brambilla <mib.mic@gmail.com>
  ///  \date Fri Nov 02 16:23:16 2012


  template <class BGF, int ORD,int DIM>
  struct StaggeredTreeLevelKernel {

    // TODO
    // ok, let's fix staggered before ..

  };
  */



}

#endif




