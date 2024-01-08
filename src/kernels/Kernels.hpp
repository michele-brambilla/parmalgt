#ifndef _KERNELS_H
#define _KERNELS_H

#include "IO.hpp"
#include "uparam.hpp"

#include <fields/LocalField.hpp>
#include <geometry/Point.hpp>
#include <background/AbelianBackground.hpp>
#include <background/TrivialBackground.hpp>
#include <background/ScalarBackground.hpp>
#include <qcd/newQCDpt.hpp>
#include <qcd/newMyQCD.hpp>
#include <types/PtTypes.hpp>

#ifdef _OPENMP
#include <omp.h>
#else
namespace kernels {
  int omp_get_max_threads() { return 1; }
  int omp_get_thread_num() { return 0; }
}
#endif

#define FLD_INFO(F) \
  typedef typename std_types<F>::ptGluon_t ptGluon;	\
  typedef typename std_types<F>::ptSU3_t ptSU3;		\
  typedef typename std_types<F>::ptsu3_t ptsu3;		\
  typedef typename std_types<F>::bgf_t BGF;		\
  typedef typename std_types<F>::point_t Point;		\
  typedef typename std_types<F>::direction_t Direction;	\
  static const int ORD = std_types<F>::order;		\
  static const int DIM = std_types<F>::n_dim;

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

namespace dirac {

  // Indices of Dirac Gamma matrices
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

  // Values of Dirac Gamma matrices
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

  // Factors multiplying each Dirac Gamma matrices to get transposed
  // ones
  const double gmT[4] = {-1.0, 1.0, -1.0, 1.0};

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
  


  template <class Field_t> 
  struct base_types {
    typedef typename Field_t::data_t data_t;
    static const int order = data_t::order;
    static const int n_dim = Field_t::dim;
    typedef pt::Point<n_dim> point_t;
    typedef pt::Direction<n_dim> direction_t;
  };



  /// struct to extract information on the field from the field type.
  template <class Field_t>
  struct std_types {
    typedef typename Field_t::data_t ptGluon_t;
    typedef typename ptGluon_t::pt_su3_t ptSU3_t;
    typedef typename ptSU3_t::pt_matrix_t ptsu3_t;
    typedef typename ptSU3_t::pt_matrix_t::SU3_t SU3_t;
    typedef typename ptGluon_t::bgf_t bgf_t;
    static const int order = ptGluon_t::order;
    static const int n_dim = Field_t::dim;
    typedef pt::Point<n_dim> point_t;
    typedef pt::Direction<n_dim> direction_t;
  };

  template <class Field_t>
  struct StapleSqKernel {
    // collect info about the field
    FLD_INFO(Field_t);

    typedef typename array_t<ptSU3, 1>::Type ptsu3_array_t;

    // checker board hyper cube size
    // c.f. geometry and localfield for more info
    static const int n_cb = 1;

    ptsu3_array_t val;
    Direction mu;

    StapleSqKernel(const Direction& nu) : mu(nu) {  }

    void operator()(const Field_t& U, const Point& n) {
      // std::cout << val[0][0] << "\n";
      val[0].zero();
      // std::cout << val[0][0] << "\n";
      for(Direction nu; nu.is_good(); ++nu)
        if(nu != mu){
          val[0] += U[n + mu][nu] *  dag(U[n][nu] * U[n + nu][mu]);
          val[0] += dag(U[n-nu][mu] * U[n+mu-nu][nu]) * U[n - nu][nu];
        }
      // Close the staple
      val[0] = U[n][mu] * val[0] ;
    }
    
    ptSU3& reduce() { 
      return val[0]; 
    }

  };

  template <class Bgf, class SUN>
  void inplace_add(SUN& s, const Bgf& b){
    for (int i = 0; i < SUN::size; ++i)
      for (int j = 0; j < SUN::size; ++j)
	s(i,j) += b(i,j);
  }
  template <class SUN>
  void inplace_add(SUN& s, const bgf::AbelianBgf& b){
    for (int i = 0; i < SUN::size; ++i)
      s(i,i) += b[i];
  }
  template <class SUN>
  void inplace_add(SUN& s, const bgf::ScalarBgf& b){
    for (int i = 0; i < SUN::size; ++i)
      s(i,i) += b.val();
  }
  // template <class SUN>
  // void inplace_add(SUN& s, const bgf::TrivialBgf&){
  //   for (int i = 0; i < SUN::size; ++i)
  //     s(i,i) += 1;
  // }
  template <class SUN, class BGF>
  void inplace_add(SUN& s, const BGF&){
    for (int i = 0; i < SUN::size; ++i)
      s(i,i) += 1;
  }
  

  ////////////////////////////////////////////////////////////
  //
  // Staple kernel taking into account the boundary correction
  // c_t. This is meant to act on the TEMPORAL components of the gauge
  // field at x_0 = 0 and x_0 = T - a.
  //
  // \author Dirk Hesse <herr.dirk.hesse@gmail.com>
  // \date Tue Oct  1 16:36:56 2013

  template <class Field_t>
  struct StapleSqKernelCtOneTemporal {
    // collect info about the field
    FLD_INFO(Field_t);

    typedef typename array_t<ptSU3, 1>::Type ptsu3_array_t;

    // checker board hyper cube size
    // c.f. geometry and localfield for more info
    static const int n_cb = 1;

    ptsu3_array_t val;
    Direction mu;

    StapleSqKernelCtOneTemporal(const Direction& nu) : mu(nu) {  }

    void operator()(Field_t& U, const Point& n) {
      static const double ct1 = -0.08896;
      val[0].zero();
      for(Direction nu; nu.is_good(); ++nu)
        if(nu != mu){
          val[0] += U[n + mu][nu] *  dag(U[n][nu] * U[n + nu][mu]);
          val[0] += dag(U[n-nu][mu] * U[n+mu-nu][nu]) * U[n - nu][nu];
        }
      // Close the staple
      val[0] = U[n][mu] * val[0];
      // Apply ct, i.e.
      // staple = staple * (1 + g_0^2 c_t^(1))
      for (int i = 2; i < ORD; ++i)
	val[0][i] += val[0][i-2] * ct1;
      // background contribution
      inplace_add(val[0][1], val[0].bgf() * ct1);
    }
    
    ptSU3& reduce() { 
      return val[0]; 
    }

  };

  ////////////////////////////////////////////////////////////
  //
  // Staple kernel taking into account the boundary correction
  // c_t. This is meant to act on the SPATIAL components of the gauge
  // field at x_0 = a.
  //
  // \author Dirk Hesse <herr.dirk.hesse@gmail.com>
  // \date Tue Oct  1 16:41:11 2013

  template <class Field_t>
  struct StapleSqKernelCtOneSpatialEarly {
    // collect info about the field
    FLD_INFO(Field_t);

    typedef typename array_t<ptSU3, 1>::Type ptsu3_array_t;

    // checker board hyper cube size
    // c.f. geometry and localfield for more info
    static const int n_cb = 1;

    ptsu3_array_t val;
    Direction mu;

    StapleSqKernelCtOneSpatialEarly(const Direction& nu) : mu(nu) {  }

    void operator()(Field_t& U, const Point& n) {
      static const double ct1 = -0.08896;
      val[0].zero();
      Direction nu;
      if(nu != mu){
	val[0] += dag(U[n-nu][mu] * U[n+mu-nu][nu]) * U[n - nu][nu];
	// contains temproal link -> factor of ct ...
	// staple = staple * (1 + g_0^2 c_t^(1))
	for (int i = 2; i < ORD; ++i)
	  val[0][i] += val[0][i-2] * ct1;
	// background contribution
	inplace_add(val[0][1], val[0].bgf() * ct1);
	val[0] += U[n + mu][nu] *  dag(U[n][nu] * U[n + nu][mu]);
      }
      for(++nu; nu.is_good(); ++nu)
        if(nu != mu){
          val[0] += U[n + mu][nu] *  dag(U[n][nu] * U[n + nu][mu]);
          val[0] += dag(U[n-nu][mu] * U[n+mu-nu][nu]) * U[n - nu][nu];
        }
      // Close the staple
      val[0] = U[n][mu] * val[0];
      // Apply ct, i.e.
    }
    
    ptSU3& reduce() { 
      return val[0]; 
    }
  };

  ////////////////////////////////////////////////////////////
  //
  // Staple kernel taking into account the boundary correction
  // c_t. This is meant to act on the SPATIAL components of the gauge
  // field at x_0 = T - a.
  //
  // \author Dirk Hesse <herr.dirk.hesse@gmail.com>
  // \date Tue Oct  1 16:42:11 2013


  template <class Field_t>
  struct StapleSqKernelCtOneSpatialLate {
    // collect info about the field
    FLD_INFO(Field_t);

    typedef typename array_t<ptSU3, 1>::Type ptsu3_array_t;

    // checker board hyper cube size
    // c.f. geometry and localfield for more info
    static const int n_cb = 1;

    ptsu3_array_t val;
    Direction mu;

    StapleSqKernelCtOneSpatialLate(const Direction& nu) : mu(nu) {  }

    void operator()(Field_t& U, const Point& n) {
      static const double ct1 = -0.08896;
      val[0].zero();
      Direction nu;
      if (nu != mu){
	val[0] += U[n + mu][nu] *  dag(U[n][nu] * U[n + nu][mu]);
	// contains temporal link -> factor of ct ...
	// staple = staple * (1 + g_0^2 c_t^(1))
	for (int i = 2; i < ORD; ++i)
	  val[0][i] += val[0][i-2] * ct1;
	// background contribution
	inplace_add(val[0][1], val[0].bgf() * ct1);
	val[0] += dag(U[n-nu][mu] * U[n+mu-nu][nu]) * U[n - nu][nu];
      }
      for(++nu; nu.is_good(); ++nu)
        if(nu != mu){
          val[0] += U[n + mu][nu] *  dag(U[n][nu] * U[n + nu][mu]);
          val[0] += dag(U[n-nu][mu] * U[n+mu-nu][nu]) * U[n - nu][nu];
        }
      // Close the staple
      val[0] = U[n][mu] * val[0];
     
    }
    
    ptSU3& reduce() { 
      return val[0]; 
    }

  };


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

  template <class Field_t>
  struct StapleReKernel {

    // collect info about the field
    typedef typename std_types<Field_t>::ptGluon_t ptGluon;
    typedef typename std_types<Field_t>::ptSU3_t ptSU3;
    typedef typename std_types<Field_t>::ptsu3_t ptsu3;
    typedef typename std_types<Field_t>::bgf_t BGF;
    typedef typename std_types<Field_t>::point_t Point;
    typedef typename std_types<Field_t>::direction_t Direction;
    static const int ORD = std_types<Field_t>::order;
    static const int DIM = std_types<Field_t>::n_dim;

    typedef typename array_t<ptSU3, 2>::Type ptsu3_array_t;    
    typedef typename array_t<double, 2>::Type weight_array_t;    
    
    // checker board hyper cube size
    // c.f. geometry and localfield for more info
    static const int n_cb = 2;

    ptsu3_array_t val;
    Direction mu;
    static weight_array_t weights;
    
    StapleReKernel(const Direction& nu) : mu(nu) {  }

    ptSU3 two_by_one(Field_t& U, const Point& n, const Direction& nu){
      return U[n + mu][mu] * U[n + mu + mu][nu] * 
             dag(U[n][nu] * U[n+nu][mu] * U[n+nu+mu][mu])
           + U[n+mu][mu] * dag( U[n-nu][mu] * 
             U[n+mu-nu][mu] * U[n+mu+mu-nu][nu] ) * U[n-nu][nu]
           + U[n+mu][nu] * dag( U[n-mu][nu] * 
             U[n-mu+nu][mu] * U[n+nu][mu] ) * U[n-mu][mu]
           + dag( U[n-nu-mu][mu] * U[n-nu][mu] * U[n-nu+mu][nu] ) * 
             U[n-nu-mu][nu] * U[n-mu][mu];
    }
    ptSU3 one_by_two(Field_t& U, const Point& n, const Direction& nu){
      return U[n + mu][nu] * U[n+mu+nu][nu] * 
             dag(U[n][nu] * U[n+nu][nu] * U[n+nu+nu][mu]) 
           + dag( U[n-nu-nu][mu] * U[n-nu-nu+mu][nu] * U[n-nu+mu][nu] ) * 
             U[n-nu-nu][nu] * U[n-nu][nu];
    }
    ptSU3 one_by_one(Field_t& U, const Point& n, const Direction& nu){
       return U[n + mu][nu] *  dag(U[n][nu] * U[n + nu][mu])
	    + dag(U[n-nu][mu] * U[n+mu-nu][nu]) * U[n - nu][nu];
    }
    
    void operator()(Field_t& U, const Point& n) {      

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

  template <class F>
  typename StapleReKernel<F>::weight_array_t StapleReKernel<F>::weights;

  template <class Field_t>
  struct TrivialPreProcess {
    static const int DIM = Field_t::dim;
    typedef typename std_types<Field_t>::point_t Point_t;
    typedef typename std_types<Field_t>::direction_t Direction_t;
    static void pre_process (Field_t& , const Point_t&, const Direction_t&) { }
    static void post_process (Field_t& , const Point_t&, const Direction_t& ) { }
  };

  // To be used at x with x_0 = a.
  // This modifies the background field, such that O(a) improvement
  // and the equations of motion hold at tree level.

  template <class Field_t>
  struct LWProcessA {
    static const int DIM = Field_t::dim;
    typedef typename std_types<Field_t>::point_t Point_t;
    typedef typename std_types<Field_t>::direction_t Direction_t;
    static void pre_process (Field_t& U, const Point_t& n, const Direction_t& k) { 
      static Direction_t t(0);
      U[n - t + k + k][t] *= 1.5;
      U[n - t - k][t] *= 1.5;
    }
    static void post_process (Field_t& U, const Point_t& n, const Direction_t& k) { 
      static Direction_t t(0);
      U[n - t + k + k][t] /= 1.5;
      U[n - t - k][t] /= 1.5;
    }
  };

  // To be used at x with x_0 = T - a.
  // This modifies the background field, such that O(a) improvement
  // and the equations of motion hold at tree level.
  template <class Field_t>
  struct LWProcessB {
    static const int DIM = Field_t::dim;
    typedef typename std_types<Field_t>::point_t Point_t;
    typedef typename std_types<Field_t>::direction_t Direction_t;
    static void pre_process (Field_t& U, const Point_t& n, const Direction_t& k) { 
      static Direction_t t(0);
      U[n - k][t] *= 1.5; // ...
      U[n + k + k][t] *= 1.5;
    }
    static void post_process (Field_t& U, const Point_t& n, const Direction_t& k) { 
      static Direction_t t(0);
      U[n - k][t] /= 1.5; // ...
      U[n + k + k][t] /= 1.5;
    }
  };


  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  ///
  ///  Kernel for the gauge update.
  ///
  ///  \tparam Field_t The type of field (I guess usually some sort of
  ///  gluon field) that the gauge update should be applied to.
  ///  \tparam StapleK_t The type of staple to use. This is employed
  ///  to implement e.g. improved gluon actions.
  ///  \tparam Process This must be a class that has two methods,
  ///  called pre_process and post_process. They are applied before
  ///  and after the gauge update is performed to be able to do things
  ///  like adjusting the weights of plaquettes at the boundary for
  ///  improved actions.
  ///
  ///  \author Dirk Hesse <herr.dirk.hesse@gmail.com>
  ///  \date Thu May 24 17:47:43 2012

  template <class Field_t, 
	    class StapleK_t, class Process >
  struct GaugeUpdateKernel {

  // collect info about the field
    typedef typename std_types<Field_t>::ptGluon_t ptGluon;
    typedef typename std_types<Field_t>::ptSU3_t ptSU3;
    typedef typename std_types<Field_t>::ptsu3_t ptsu3;
    typedef typename std_types<Field_t>::bgf_t BGF;
    typedef typename std_types<Field_t>::point_t Point;
    typedef typename std_types<Field_t>::direction_t Direction;
    static const int ORD = std_types<Field_t>::order;
    static const int DIM = std_types<Field_t>::n_dim;

    typedef std::vector<Cplx>::iterator cpx_vec_it;
    typedef std::vector<std::vector<Cplx> >::iterator outer_cvec_it;

    // checker board hyper cube size
    // c.f. geometry and localfield for more info
    static const int n_cb = StapleK_t::n_cb;    

    // for testing, c.f. below
    static std::vector<ranlxd::Rand> rands;
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



    void operator()(Field_t& U, const Point& n) {
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
      //
      // WARNING
      // since we use stau as a pre-factor, it looks like the random
      // matrices are normalized to sqrt(2)!
      tmp[0] -= stau*sun::SU3rand(rands.at(n));
      U[n][mu] = exp<BGF, ORD>(tmp)*U[n][mu]; // back to SU3
      //#pragma omp critical // TODO maybe one can use a reduce or so here
#ifndef SF
      M[omp_get_thread_num()] += get_q(U[n][mu]); // zero momentum contribution
#endif
    }
    
    void reduce(){
      //outer_cvec_it i = plaq.begin();
      //outer_cvec_it pe = plaq.end();
      //if (++i != plaq.end())
      //  for (; i != pe; ++i)
      //	  plaq[0] += *i;
#ifndef SF
      typename std::vector<ptsu3>::iterator j = M.begin();
      if (++j != M.end())
        for (; j != M.end(); ++j)
          M[0] += *j;
#endif
    }
    
  };
  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  ///
  ///  Kernel for the Wilson flow, simple Euler scheme.
  ///
  ///  \tparam Field_t The type of field (I guess usually some sort of
  ///  gluon field) that the gauge update should be applied to.
  ///  \tparam StapleK_t The type of staple to use. This is employed
  ///  to implement e.g. improved gluon actions.
  ///  \tparam Process This must be a class that has two methods,
  ///  called pre_process and post_process. They are applied before
  ///  and after the gauge update is performed to be able to do things
  ///  like adjusting the weights of plaquettes at the boundary for
  ///  improved actions.
  ///
  ///  \author Dirk Hesse <herr.dirk.hesse@gmail.com>
  ///  \date Thu May 24 17:47:43 2012

  template <class Field_t, 
	    class StapleK_t, class Process >
  struct WilFlowKernel {

  // collect info about the field
    typedef typename std_types<Field_t>::ptGluon_t ptGluon;
    typedef typename std_types<Field_t>::ptSU3_t ptSU3;
    typedef typename std_types<Field_t>::ptsu3_t ptsu3;
    typedef typename std_types<Field_t>::bgf_t BGF;
    typedef typename std_types<Field_t>::point_t Point;
    typedef typename std_types<Field_t>::direction_t Direction;
    static const int ORD = std_types<Field_t>::order;
    static const int DIM = std_types<Field_t>::n_dim;

    // checker board hyper cube size
    // c.f. geometry and localfield for more info
    static const int n_cb = StapleK_t::n_cb;    

    Direction mu;

    double taug;

    WilFlowKernel(const Direction& nu, const double& t) :
      mu(nu), taug(t) { }

    void operator()(Field_t& U, const Point& n) {
      // Make a Kernel to calculate and store the plaquette(s)
      StapleK_t st(mu); // maye make a vector of this a class member
      Process::pre_process(U,n,mu);
      st(U,n);
      Process::post_process(U,n,mu);
      
      ptsu3 tmp  = st.reduce().reH() * -taug;

      U[n][mu] = exp<BGF, ORD>(tmp)*U[n][mu]; // back to SU3

    }
  };
  


  template <class Field_t, 
	    class StapleK_t, class Process >
  struct WilFlowMeasKernel {

  // collect info about the field
    typedef typename std_types<Field_t>::ptGluon_t ptGluon;
    typedef typename std_types<Field_t>::ptSU3_t ptSU3;
    typedef typename std_types<Field_t>::ptsu3_t ptsu3;
    typedef typename std_types<Field_t>::bgf_t BGF;
    typedef typename std_types<Field_t>::point_t Point;
    typedef typename std_types<Field_t>::direction_t Direction;
    static const int ORD = std_types<Field_t>::order;
    static const int DIM = std_types<Field_t>::n_dim;

    // checker board hyper cube size
    // c.f. geometry and localfield for more info
    static const int n_cb = StapleK_t::n_cb;    

    Direction mu;

    double taug, stau;

    WilFlowMeasKernel(const Direction& nu, const double& t,
                  Field_t& targ) :
      mu(nu), taug(t), stau(sqrt(t)), target(&targ) { }



    void operator()(Field_t& U, const Point& n) {
      ptSU3 W;

      // Make a Kernel to calculate and store the plaquette(s)
      StapleK_t st(mu); // maye make a vector of this a class member
      Process::pre_process(U,n,mu);
      st(U,n);
      Process::post_process(U,n,mu);
      
      ptsu3 tmp  = st.reduce().reH() * -taug;

      (*target)[n][mu] = exp<BGF, ORD>(tmp); // back to SU3

    }

    Field_t *target;
  };

  template <class Field_t,
	    class StapleK_t, class Process >
  struct WilFlowApplyKernel {

  // collect info about the field
    typedef typename std_types<Field_t>::ptGluon_t ptGluon;
    typedef typename std_types<Field_t>::ptSU3_t ptSU3;
    typedef typename std_types<Field_t>::ptsu3_t ptsu3;
    typedef typename std_types<Field_t>::bgf_t BGF;
    typedef typename std_types<Field_t>::point_t Point;
    typedef typename std_types<Field_t>::direction_t Direction;
    static const int ORD = std_types<Field_t>::order;
    static const int DIM = std_types<Field_t>::n_dim;

    // checker board hyper cube size
    // c.f. geometry and localfield for more info
    static const int n_cb = StapleK_t::n_cb;

    Direction mu;

    double taug, stau;

    WilFlowApplyKernel(const Direction& nu, const double& t,
                       Field_t &src) :
      mu(nu), taug(t), stau(sqrt(t)), source(&src) { }


    void operator()(Field_t& U, const Point& n) {
      ptSU3 W;

      // Make a Kernel to calculate and store the plaquette(s)
      StapleK_t st(mu); // maye make a vector of this a class member
      Process::pre_process(U,n,mu);
      st(U,n);
      Process::post_process(U,n,mu);

      ptsu3 tmp  = st.reduce().reH() * -taug;

      U[n][mu] = (*source)[n][mu]*U[n][mu]; // back to SU3
    }
    Field_t *source;
  };


template <class C, class P, class Q> std::vector<ranlxd::Rand>
  kernels::GaugeUpdateKernel<C,P,Q>::rands;

  template <class Field_t> struct RSU3Kernel {

    // collect info about the field
    static const int n_dim = Field_t::dim;
    typedef pt::Point<n_dim> Point;
    // Choose random number generator to use
    // Note that you have to change the same lines
    // in Methods.hpp
    //*/
    // using RANLUX
    using rand_vec_t = std::vector<ranlxd::Rand>;
    /*/
    // using MyRand
    typedef std::vector<MyRand> rand_vec_t;
    //*/
    static rand_vec_t rands;
    static const int n_cb = 0;

    void operator()(Field_t& U, const Point& n) const {
      U[n] = sun::SU3rand(rands[n]);
    }
  };

  template <class C> typename kernels::RSU3Kernel<C>::rand_vec_t
  kernels::RSU3Kernel<C>::rands;


  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  ///
  ///  Kernel for gauge fixing.
  ///
  ///  Note I implemented three different methods of which only the
  ///  third one seems to work fine at the moment.
  ///
  ///  \tparam METHOD  Gauge fixing mode. WARNING: Only mode 1 and 3
  ///  are tested and only 1 is fully trusted!
  ///  \tparam Field_t The kind of field to be used.
  ///
  ///  \author Dirk Hesse <herr.dirk.hesse@gmail.com>
  ///  \date Thu May 24 17:49:09 2012
  ///

  template <int METHOD, class Field_t>
  class GaugeFixingKernel {
  public:
    // collect info about the field
    typedef typename std_types<Field_t>::ptGluon_t ptGluon;
    typedef typename std_types<Field_t>::ptSU3_t ptSU3;
    typedef typename std_types<Field_t>::ptsu3_t ptsu3;
    typedef typename std_types<Field_t>::bgf_t BGF;
    typedef typename std_types<Field_t>::point_t Point;
    typedef typename std_types<Field_t>::direction_t Direction;
    static const int ORD = std_types<Field_t>::order;
    static const int DIM = std_types<Field_t>::n_dim;

    explicit GaugeFixingKernel (const double& a) : alpha (a) { }
    void operator()(Field_t& U, const Point& n) const { 
    do_it(U, n, mode_selektor<METHOD>());
    }

    // checker board hyper cube size
    // c.f. geometry and localfield for more info
    static const int n_cb = 1;

private:
    double alpha;
  template <int M> struct mode_selektor { };
  void do_it(Field_t& U, const Point& n, 
             const mode_selektor<1>&) const {
    // exp version
    ptSU3 omega;
    for (Direction mu; mu.is_good(); ++mu)
      omega += U[n - mu][mu] - U[n][mu] ;
    
    ptSU3 Omega = exp<BGF, ORD>( alpha * omega.reH());
    ptSU3 OmegaDag = exp<BGF, ORD>( -alpha * omega.reH());
    
    for (Direction mu; mu.is_good(); ++mu){
      U[n][mu] = Omega * U[n][mu];
      U[n - mu][mu] *= OmegaDag;
    }
  }
  void do_it(Field_t& U, const Point& n,
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
  void do_it(Field_t& U, const Point& n,
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

  void do_it(Field_t& U, const Point& n, 
             const mode_selektor<4>&) const {
    // exp version
    ptSU3 omega;
    for (Direction mu; mu.is_good(); ++mu)
      omega += dag(U[n-mu][mu].bgf())*U[n - mu][mu] 
	- U[n][mu]*dag(U[n][mu].bgf());
    
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

  template <class Field_t>
  struct ZeroModeSubtractionKernel
  {

    // collect info about the field
    typedef typename std_types<Field_t>::ptGluon_t ptGluon;
    typedef typename std_types<Field_t>::ptSU3_t ptSU3;
    typedef typename std_types<Field_t>::ptsu3_t ptsu3;
    typedef typename std_types<Field_t>::bgf_t BGF;
    typedef typename std_types<Field_t>::point_t Point;
    typedef typename std_types<Field_t>::direction_t Direction;
    static const int ORD = std_types<Field_t>::order;
    static const int DIM = std_types<Field_t>::n_dim;

    // checker board hyper cube size
    // c.f. geometry and localfield for more info
    static const int n_cb = 1;

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

    void operator()(Field_t& U, const Point& n) {
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
  template <class Field_t>
  struct SetBgfKernel {

    // collect info about the field
    typedef typename std_types<Field_t>::ptGluon_t ptGluon;
    typedef typename std_types<Field_t>::ptSU3_t ptSU3;
    typedef typename std_types<Field_t>::ptsu3_t ptsu3;
    typedef typename std_types<Field_t>::bgf_t BGF;
    typedef typename std_types<Field_t>::point_t Point;
    typedef typename std_types<Field_t>::direction_t Direction;
    static const int ORD = std_types<Field_t>::order;
    static const int DIM = std_types<Field_t>::n_dim;

    // checker board hyper cube size
    // c.f. geometry and localfield for more info
    static const int n_cb = 0;

    explicit SetBgfKernel(const int& t_in) : t(t_in) { }
    int t;
    void operator()(Field_t& U, const Point& n) const {
      impl(U, n, BGF());
    }

    void impl(Field_t& U, const Point& n, const bgf::AbelianBgf&) const {
      for (Direction mu; mu.is_good(); ++mu)
	//#ifdef HIGHER_ORDER_INT
        //U[n][mu].bgf() =  bgf::get_abelian_bgf(t, 0);
	//#else
        U[n][mu].bgf() =  bgf::get_abelian_bgf(t, mu);
      //#endif
    }

    void impl(Field_t& U, const Point& n, const bgf::ScalarBgf&)
    const { }

  };


  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  ///
  ///  Kernel to measure the norm.
  ///
  ///  \author Dirk Hesse <herr.dirk.hesse@gmail.com>
  ///  \date Thu May 24 17:52:56 2012
  template <class Field_t>
  struct MeasureNormKernel {

    // collect info about the field
    typedef typename std_types<Field_t>::ptGluon_t ptGluon;
    typedef typename std_types<Field_t>::ptSU3_t ptSU3;
    typedef typename std_types<Field_t>::ptsu3_t ptsu3;
    typedef typename std_types<Field_t>::bgf_t BGF;
    typedef typename std_types<Field_t>::point_t Point;
    typedef typename std_types<Field_t>::direction_t Direction;
    static const int ORD = std_types<Field_t>::order;
    static const int DIM = std_types<Field_t>::n_dim;

    // checker board hyper cube size
    // c.f. geometry and localfield for more info
    static const int n_cb = 0;
    std::vector<typename array_t<double, ORD+1>::Type> norm;

    explicit MeasureNormKernel() : norm(omp_get_max_threads()) { }

    void operator()(const Field_t& U, const Point& n) {
      std::vector<double> tmp = U[n].Norm();
      int i = omp_get_thread_num();
      for (int k = 0; k < ORD+1; ++k)
        norm[i][k] += tmp[k];
    }
    typename array_t<double, ORD+1>::Type reduce(){
      for (int i = 1, j = omp_get_max_threads(); i < j; ++i)
        for (int k = 0; k < ORD+1; ++k)
          norm[0][k] += norm[i][k];
      return norm[0];
    }
    typename array_t<double, ORD+1>::Type reduce(typename array_t<double, ORD+1>::Type& other){
      other = norm[0];
      for (int i = 1, j = omp_get_max_threads(); i < j; ++i)
        for (int k = 0; k < ORD+1; ++k) {
          norm[0][k] += norm[i][k];
	  other[k]   += norm[i][k];
	}
      return norm[0];
    }
    typename array_t<double, ORD+1>::Type reduce(Norm<ORD+1>& other){
      for (int k = 0; k < ORD+1; ++k)
	other[k] = norm[0][k];
      for (int i = 1, j = omp_get_max_threads(); i < j; ++i)
        for (int k = 0; k < ORD+1; ++k) {
          norm[0][k] += norm[i][k];
	  other[k] += norm[i][k];
	}
      return norm[0];
    }
  };

  // measures Udagger * U
  template <class Field_t>
  struct UdagUKernel {

    // collect info about the field
    typedef typename std_types<Field_t>::ptGluon_t ptGluon;
    typedef typename std_types<Field_t>::ptSU3_t ptSU3;
    typedef typename std_types<Field_t>::ptsu3_t ptsu3;
    typedef typename std_types<Field_t>::bgf_t BGF;
    typedef typename std_types<Field_t>::point_t Point;
    typedef typename std_types<Field_t>::direction_t Direction;
    static const int ORD = std_types<Field_t>::order;
    static const int DIM = std_types<Field_t>::n_dim;

    // checker board hyper cube size
    // c.f. geometry and localfield for more info
    static const int n_cb = 0;

    ptSU3 val;
    UdagUKernel () : val(bgf::zero<BGF>()) { }
    void operator()(Field_t& U, const Point& n){
      ptSU3 tmp(bgf::zero<BGF>());
      for (Direction mu; mu.is_good(); ++mu)
        tmp += dag(U[n][mu]) * U[n][mu];
#pragma omp critical
      val += tmp;
    }
  
};

  template <class Field_t, class OutputIterator>
  struct Buffer {
    // collect info about the field
    static const int n_dim = Field_t::dim;
    typedef pt::Point<n_dim> point_t;
    typedef pt::Direction<n_dim> direction_t;

    // This may NOT be executed in parallel, so ...
    typedef void NoPar;
    
    // iterator to the positon where the next write goes
    OutputIterator oi;
    
    Buffer(OutputIterator o) : oi(o) { }

    void operator()(Field_t& U, const point_t& n) {
      *(oi++) = U[n];
    }
  };
  template <class Field_t, class InputIterator>
  struct Unbuffer {
    // collect info about the field
    static const int n_dim = Field_t::dim;
    typedef pt::Point<n_dim> point_t;
    typedef pt::Direction<n_dim> direction_t;

    // This may NOT be executed in parallel, so ...
    typedef void NoPar;
    
    // iterator to the positon where the next write goes
    InputIterator ii;
    
    Unbuffer(InputIterator i) : ii(i) { }

    void operator()(Field_t& U, const point_t& n){
      U[n] = *(ii++);
    }
  };

// Measure the average paquette
  template <class Field_t>
  struct PlaqKernel {

    // collect info about the field
    typedef typename std_types<Field_t>::ptGluon_t ptGluon;
    typedef typename std_types<Field_t>::ptSU3_t ptSU3;
    typedef typename std_types<Field_t>::ptsu3_t ptsu3;
    typedef typename std_types<Field_t>::bgf_t BGF;
    typedef typename std_types<Field_t>::point_t Point;
    typedef typename std_types<Field_t>::direction_t Direction;
    static const int ORD = std_types<Field_t>::order;
    static const int DIM = std_types<Field_t>::n_dim;

    // checker board hyper cube size
    // c.f. geometry and localfield for more info
    static const int n_cb = 0;

    ptSU3 val;
    
    PlaqKernel () : val(bgf::zero<BGF>()) { }
    
    void operator()(Field_t& U, const Point& n){
      ptSU3 tmp(bgf::zero<BGF>());
      for (Direction mu(0); mu.is_good(); ++mu)
        for (Direction nu(0); nu.is_good(); ++nu)
          tmp += U[n][mu] * U[n + mu][nu] * 
            dag( U[n][nu] * U[n + nu][mu] );
#pragma omp critical
      val += tmp;
    }
  };
  
  // Kernel to construct the gauge fixing function at t = 0
  template <class Field_t>
  struct GFMeasKernel {

    // collect info about the field
    typedef typename std_types<Field_t>::ptGluon_t ptGluon;
    typedef typename std_types<Field_t>::ptSU3_t ptSU3;
    typedef typename std_types<Field_t>::ptsu3_t ptsu3;
    typedef typename std_types<Field_t>::bgf_t BGF;
    typedef typename std_types<Field_t>::point_t Point;
    typedef typename std_types<Field_t>::direction_t Direction;
    static const int ORD = std_types<Field_t>::order;
    static const int DIM = std_types<Field_t>::n_dim;

    // checker board hyper cube size
    // c.f. geometry and localfield for more info
    static const int n_cb = 0;

    ptsu3 val;

    void operator()(Field_t& U, const Point& n){
#pragma omp critical
      val += get_q(U[n][Direction(0)]);
    }
    
  };  
  // Kernel to execute the gauge fixing function at t = 0
  // (the one where we used U_0(\vec y, 0) to construct the gf 
  // function)
  template <class Field_t>
  struct GFApplyKernel {

    // collect info about the field
    typedef typename std_types<Field_t>::ptGluon_t ptGluon;
    typedef typename std_types<Field_t>::ptSU3_t ptSU3;
    typedef typename std_types<Field_t>::ptsu3_t ptsu3;
    typedef typename std_types<Field_t>::bgf_t BGF;
    typedef typename std_types<Field_t>::point_t Point;
    typedef typename std_types<Field_t>::direction_t Direction;
    static const int ORD = std_types<Field_t>::order;
    static const int DIM = std_types<Field_t>::n_dim;

    ptSU3 Omega, OmegaDag;
    
    // checker board hyper cube size
    // c.f. geometry and localfield for more info
    static const int n_cb = 0;

    GFApplyKernel (ptsu3 omega, const double& alpha,
                                 const int& L) { 
      for (int r = 0; r < ORD; ++r)
        for (int i = 0; i < 3; ++i)
          for (int j = 0; j < 3; ++j)
            if (i != j)
              omega[r](i,j) = 0;
      
      omega /= L*L*L;
      Omega = exp<BGF, ORD>( -alpha * omega.reH());
    }

    void operator()(Field_t& U, const Point& n){
      static Direction t(0);
      U[n][t] = Omega * U[n][t];
    }
    
  };
  template <class Field_t>
  struct GFApplyKernelTRBG {

    // collect info about the field
    typedef typename std_types<Field_t>::ptGluon_t ptGluon;
    typedef typename std_types<Field_t>::ptSU3_t ptSU3;
    typedef typename std_types<Field_t>::ptsu3_t ptsu3;
    typedef typename std_types<Field_t>::bgf_t BGF;
    typedef typename std_types<Field_t>::point_t Point;
    typedef typename std_types<Field_t>::direction_t Direction;
    static const int ORD = std_types<Field_t>::order;
    static const int DIM = std_types<Field_t>::n_dim;

    ptSU3 Omega, OmegaDag;
    
    // checker board hyper cube size
    // c.f. geometry and localfield for more info
    static const int n_cb = 0;

    GFApplyKernelTRBG (ptsu3 omega, const double& alpha,
                                 const int& L) { 
      omega /= L*L*L;
      Omega = exp<BGF, ORD>( -alpha * omega.reH());
    }

    void operator()(Field_t& U, const Point& n){
      static Direction t(0);
      U[n][t] = Omega * U[n][t];
    }
    
  };

  // helper class to get around compile errors during initialization
  // of the kernels below ..
  template <class B> struct init_helper_gamma {
    void operator()(const int&, B&) const {}
  };
  template <> 
  struct init_helper_gamma<bgf::AbelianBgf> {
    void operator()(const int& L, bgf::AbelianBgf& Ctilde) const {
      Cplx ioL(0, 1./L);
      Ctilde[0] = -2.*ioL;
      Ctilde[1] = ioL;
      Ctilde[2] = ioL;
    }
  };
  template <class B> struct init_helper_vbar {
    void operator()(const int&, B&) const {}
  };
  template <> 
  struct init_helper_vbar<bgf::AbelianBgf> {
    void operator()(const int& L, bgf::AbelianBgf& Ctilde) const {
      Cplx ioL(0, 1./L);
      Ctilde[0] = 0;
      Ctilde[1] = 2.*ioL;
      Ctilde[2] = -2.*ioL;
    }
  };

#ifdef IMP_ACT
  // Measure Gamma at t = 0
  template <class Field_t, template <class C> class init_helper_t>
  struct GammaLowerKernel {

    // collect info about the field
    typedef typename std_types<Field_t>::ptGluon_t ptGluon;
    typedef typename std_types<Field_t>::ptSU3_t ptSU3;
    typedef typename std_types<Field_t>::ptsu3_t ptsu3;
    typedef typename std_types<Field_t>::bgf_t BGF;
    typedef typename std_types<Field_t>::point_t Point;
    typedef typename std_types<Field_t>::direction_t Direction;
    static const int ORD = std_types<Field_t>::order;
    static const int DIM = std_types<Field_t>::n_dim;

    // checker board hyper cube size
    // c.f. geometry and localfield for more info
    static const int n_cb = 0;

    // improvement coefficents
    typedef typename array_t<double, 2>::Type weights;
    weights c;

    ptSU3 val;
    // \tilde C = - [d_eta C]
    // at the t = 0 side, we have dagger(e^C) and hence an insertion
    // of \tilde C, since C itself is purely imaginary
    BGF Ctilde;
    explicit GammaLowerKernel (const int& L, const weights& c_) : val(bgf::zero<BGF>()), c(c_) { 
      init_helper_t<BGF>()(L, Ctilde);
    }
    void operator()(const Field_t& U, const Point& n){
      Direction t(0);
      ptSU3 tmp(bgf::zero<BGF>());
      for (Direction k(1); k.is_good(); ++k){
        // the 1x1 contribution
        tmp += (Ctilde * dag(U[n][k]) * U[n][t] * 
                U[n + t][k] * dag( U[n + k][t] )) * c[0];
        // the 2x1 contribution w\ two links @ the boundary
        // we actually have two contributions to the derivative form
        // each boundary link. However, we do not care because the
        // insertion of \tilde C commutes with the links. Hence, we
        // just multiply the weight with two, 
        //  w = 2 * (3/2) * c_1 = 3*(-1/12) = -1/4
        tmp += (Ctilde * dag(U[n][k]) * dag(U[n-k][k]) 
                * U[n-k][t] * U[n + t - k][k]
                * U[n + t][k] * dag( U[n + k][t] )) * 2. * (3./2.) * c[1];
        // the 1x2 contribution, with usual weight c_1
        tmp += (Ctilde * dag(U[n][k]) * U[n][t] 
                * U[n+t][t] * U[n+t+t][k] 
                * dag(U[n + k + t][t]) * dag(U[n + k][t])) * c[1];
      }
#pragma omp critical
      val += tmp;
    }
  };

  // Measure Gamma at t = T - a
  template <class Field_t, template <class C> class init_helper_t>
  struct GammaUpperKernel {

    // collect info about the field
    typedef typename std_types<Field_t>::ptGluon_t ptGluon;
    typedef typename std_types<Field_t>::ptSU3_t ptSU3;
    typedef typename std_types<Field_t>::ptsu3_t ptsu3;
    typedef typename std_types<Field_t>::bgf_t BGF;
    typedef typename std_types<Field_t>::point_t Point;
    typedef typename std_types<Field_t>::direction_t Direction;
    static const int ORD = std_types<Field_t>::order;
    static const int DIM = std_types<Field_t>::n_dim;

    // checker board hyper cube size
    // c.f. geometry and localfield for more info
    static const int n_cb = 0;

    // improvement coefficents
    typedef typename array_t<double, 2>::Type weights;
    weights c;

    ptSU3 val;
    // here, we need [d_eta C'], which is equal to -[d_eta C], hence
    // we can use Ctilde as above
    BGF Ctilde;
    explicit GammaUpperKernel (const int& L, const weights& c_) : val(bgf::zero<BGF>()), c(c_) { 
      init_helper_t<BGF>()(L, Ctilde);
    }
    void operator()(const Field_t& U, const Point& n){
      Direction t(0);
      ptSU3 tmp(bgf::zero<BGF>());
      for (Direction k(1); k.is_good(); ++k){
        // the 1x1 contribution
        tmp += (Ctilde * U[n + t][k] * dag(U[n + k][t])
                * dag(U[n][k]) * U[n][t]) * c[0];
        // the 2x1 contribution w\ two links @ the boundary
        // we actually have two contributions to the derivative form
        // each boundary link. However, we do not care because the
        // insertion of \tilde C commutes with the links. Hence, we
        // just multiply the weight with two, 
        //  w = 2 * (3/2) * c_1 = 3*(-1/12) = -1/4
        tmp += (Ctilde * U[n + t][k] * U[n + t + k][k]
                * dag(U[n + k + k][t]) * dag(U[n + k][k])
                * dag(U[n][k]) * U[n][t]) * 2. * (3./2.) * c[1];
        // the 1x2 contribution, with usual weight c_1
        tmp += (Ctilde * U[n + t][k] * dag(U[n + k][t])
                * dag(U[n + k - t][t]) * dag(U[n-t][k]) 
                * U[n - t][t] * U[n][t]) * c[1];
      }
#pragma omp critical
      val += tmp;
    }
  };

#else
  // Measure Gamma at t = 0
  template <class Field_t, template <class C> class init_helper_t>
  struct GammaLowerKernel {

    // collect info about the field
    typedef typename std_types<Field_t>::ptGluon_t ptGluon;
    typedef typename std_types<Field_t>::ptSU3_t ptSU3;
    typedef typename std_types<Field_t>::ptsu3_t ptsu3;
    typedef typename std_types<Field_t>::bgf_t BGF;
    typedef typename std_types<Field_t>::point_t Point;
    typedef typename std_types<Field_t>::direction_t Direction;
    static const int ORD = std_types<Field_t>::order;
    static const int DIM = std_types<Field_t>::n_dim;

    // checker board hyper cube size
    // c.f. geometry and localfield for more info
    static const int n_cb = 0;

    ptSU3 val;
    // \tilde C = - [d_eta C]
    // at the t = 0 side, we have dagger(e^C) and hence an insertion
    // of \tilde C, since C itself is purely imaginary
    BGF Ctilde;
    explicit GammaLowerKernel (const int &L) : val(bgf::zero<BGF>()) { 
      init_helper_t<BGF>()(L, Ctilde);
    }
    void operator()(const Field_t& U, const Point& n){
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
  template <class Field_t, template <class C> class init_helper_t>
  struct GammaUpperKernel {

    // collect info about the field
    typedef typename std_types<Field_t>::ptGluon_t ptGluon;
    typedef typename std_types<Field_t>::ptSU3_t ptSU3;
    typedef typename std_types<Field_t>::ptsu3_t ptsu3;
    typedef typename std_types<Field_t>::bgf_t BGF;
    typedef typename std_types<Field_t>::point_t Point;
    typedef typename std_types<Field_t>::direction_t Direction;
    static const int ORD = std_types<Field_t>::order;
    static const int DIM = std_types<Field_t>::n_dim;

    // checker board hyper cube size
    // c.f. geometry and localfield for more info
    static const int n_cb = 0;

    ptSU3 val;
    // here, we need [d_eta C'], which is equal to -[d_eta C], hence
    // we can use Ctilde as above
    BGF Ctilde;
    explicit GammaUpperKernel (const int& L) : val(bgf::zero<BGF>()) { 
      init_helper_t<BGF>()(L, Ctilde);
    }
    void operator()(const Field_t& U, const Point& n){
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


  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  ///
  ///  Writing a gluon to a file.
  ///
  ///  \author Dirk Hesse <herr.dirk.hesse@gmail.com>
  ///  \date Fri May 25 15:59:06 2012

  template <class Field_t>
  class FileWriterKernel {
  public:
    // collect info about the field
    typedef typename std_types<Field_t>::ptGluon_t ptGluon;
    typedef typename std_types<Field_t>::ptSU3_t ptSU3;
    typedef typename std_types<Field_t>::ptsu3_t ptsu3;
    typedef typename std_types<Field_t>::bgf_t BGF;
    typedef typename std_types<Field_t>::point_t Point;
    typedef typename std_types<Field_t>::direction_t Direction;
    static const int ORD = std_types<Field_t>::order;
    static const int DIM = std_types<Field_t>::n_dim;

    // This may NOT be executed in parallel, so ...
    typedef void NoPar;

    explicit FileWriterKernel (uparam::Param& p) : o(p) { }

    void operator()(Field_t& U, const Point& n){
      for (Direction mu(0); mu.is_good(); ++mu)
#pragma omp critical
        U[n][mu].write(o);
    }
    io::CheckedOut o;
  private:
    // make n_cb private to prevent parallel application of this
    // kernel, because this would be a terrible idea
    static const int n_cb = 0;
  };
  
  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  ///
  ///  Reading a gluon from a file.
  ///
  ///  \author Dirk Hesse <herr.dirk.hesse@gmail.com>
  ///  \date Wed May 30 18:37:03 2012

  template <class Field_t>
  struct FileReaderKernel {
  public:
    // collect info about the field
    typedef typename std_types<Field_t>::ptGluon_t ptGluon;
    typedef typename std_types<Field_t>::ptSU3_t ptSU3;
    typedef typename std_types<Field_t>::ptsu3_t ptsu3;
    typedef typename std_types<Field_t>::bgf_t BGF;
    typedef typename std_types<Field_t>::point_t Point;
    typedef typename std_types<Field_t>::direction_t Direction;
    static const int ORD = std_types<Field_t>::order;
    static const int DIM = std_types<Field_t>::n_dim;

    // This may NOT be executed in parallel, so ...
    typedef void NoPar;

    explicit FileReaderKernel (uparam::Param& p) : i(p) { }

    void operator()(Field_t& U, const Point& n){
      for (Direction mu(0); mu.is_good(); ++mu)
        U[n][mu].read(i);
    }
    io::CheckedIn i;
  private:
    // make n_cb private to prevent parallel application of this
    // kernel, because this would be a terrible idea
    static const int n_cb = 0;
  };




  namespace detail {
    
    template<class Field_t>
    struct Mass : base_types<Field_t> {
      typedef typename base_types<Field_t>::data_t F;
      Mass() { }
      void operator()(const F& src, F& dest,const double& mass) const {
	dest = src*mass;
      }
    };
    template<class Field_t>
    struct Mass5 : base_types<Field_t> {
      typedef typename base_types<Field_t>::data_t F;
      typedef typename base_types<Field_t>::direction_t Direction;
      Mass5() { }
      void operator()(const F& src, F& dest,const double& mass) const {
	for(Direction nu(0);nu.is_good();++nu)
	  dest[nu] = (src[dirac::gmuind[4][nu]]*dirac::gmuval[4][nu]*mass);
      }
    };


    template<class Gauge_t, class Field_t>
    struct Gamma : public base_types<Field_t> {
      typedef typename base_types<Field_t>::data_t F;
      typedef typename base_types<Field_t>::direction_t Direction;
      typedef typename base_types<Field_t>::point_t Point;
      Gamma() { }
      void operator()(const Field_t& src,const Point& n,const Direction& mu) {
	Point dn = n-Direction(mu);
	Point up = n+Direction(mu);
	for(Direction nu(0);nu.is_good();++nu) {
	  Xi1[nu] = (src[dn][nu]+src[dn][dirac::gmuind[mu][nu]]*dirac::gmuval[mu][nu]);
	  Xi2[nu] = (src[up][nu]-src[up][dirac::gmuind[mu][nu]]*dirac::gmuval[mu][nu]);
	}
      }
      const F& xi1() { return Xi1; }
      const F& xi2() { return Xi2; }
    private:
      F Xi1, Xi2;
    };
    template<class Gauge_t, class Field_t>
    struct Gamma5 : public base_types<Field_t> {
      typedef typename base_types<Field_t>::data_t F;
      typedef typename base_types<Field_t>::direction_t Direction;
      typedef typename base_types<Field_t>::point_t Point;
      Gamma5() { }
      void operator()(const Field_t& src,const Point& n,const Direction& mu) {
	Point dn = n-Direction(mu);
	Point up = n+Direction(mu);
	for(Direction nu(0);nu.is_good();++nu) {
	  Xi1[nu] = (src[dn][dirac::gmuind[4        ][nu]]*dirac::gmuval[4        ][nu] + 
		     src[dn][dirac::gmuind[5+(int)mu][nu]]*dirac::gmuval[5+(int)mu][nu] );
	  Xi2[nu] = (src[up][dirac::gmuind[4        ][nu]]*dirac::gmuval[4        ][nu] - 
		     src[up][dirac::gmuind[5+(int)mu][nu]]*dirac::gmuval[5+(int)mu][nu] );
	}
      }
      const F& xi1() { return Xi1; }
      const F& xi2() { return Xi2; }
    private:
      F Xi1, Xi2;
    };
    
    
      

  }


  enum boundaryId{bulk,lower,upper};

  // template <class Field_t, int boundary, template<class> class Mass_t, template<class,class> class Gamma_t>
  // struct WilsonTLBase : public base_types<Field_t> {
  // public:
  //   // collect info about the field
  //   typedef typename std_types<Field_t>::point_t Point;
  //   typedef typename std_types<Field_t>::direction_t Direction;
  //   static const int DIM = std_types<Field_t>::n_dim;
        
  //   // fermion
  //   typedef SpinColor<DIM> Fermion;
  //   typedef typename fields::LocalField<Fermion, DIM> FermionField;

  //   // checker board hyper cube size
  //   // c.f. geometry and localfield for more info
  //   static const int n_cb = 0;
    
  //   WilsonTLBase(const Field_t& G, const FermionField& X, 
  // 			   const double& m ) :
  //     U(G), src(X), mass(m) { };
    
  //   void operator() ( FermionField& dest, const Point& n) {
  //     M(src[n],dest[n],mass);
  //     do_spatial(dest,n);
  //     do_temporal(dest,n,mode_selektor<boundary>());
  //   }
  
  // private:
    
  //   const double& mass;
  //   const Field_t& U;
  //   const FermionField& src;
  //   Mass_t<FermionField> M;
  //   Gamma_t<Field_t,FermionField> gmm;

  //   template <int M> struct mode_selektor { };
  
  //   void do_spatial(FermionField& dest,const Point& n) {
  //     for( Direction mu(1); mu.is_good(); ++mu ) {
  // 	gmm(src,n,mu);
  // 	dest[n] -= (dag(U[n-mu][mu].bgf())*gmm.xi1() +
  // 		    U[n ][mu].bgf()*gmm.xi2())*.5;
  //     }
  //   }
  //   void do_temporal(FermionField& dest,const Point& n,
  //   		     const mode_selektor<bulk>) {
  //     gmm(src,n,Direction(0));
  //     dest[n] -= (dag(U[n-Direction(0)][Direction(0)].bgf())*gmm.xi1() +
  //   		  U[n ][Direction(0)].bgf()*gmm.xi2())*.5;
  //   }
  //   void do_temporal(FermionField& dest,const Point& n,
  //   		     const mode_selektor<lower>) {
  //     gmm(src,n,Direction(0));
  //     dest[n] += (dag(U[n-Direction(0)][Direction(0)].bgf())*gmm.xi1() -
  //   		  U[n ][Direction(0)].bgf()*gmm.xi2())*.5;
  //   }
  //   void do_temporal(FermionField& dest,const Point& n,
  //   		     const mode_selektor<upper>) {
  //     gmm(src,n,Direction(0));
  //     dest[n] -= (dag(U[n-Direction(0)][Direction(0)].bgf())*gmm.xi1() -
  //   		  U[n ][Direction(0)].bgf()*gmm.xi2())*.5;
  //   }  

  // };

  // template <class Field_t, int boundary, template<class> class Mass_t, template<class,class> class Gamma_t>
  // struct WilsonBase : public base_types<Field_t> {
  // public:
  //   // collect info about the field
  //   typedef typename std_types<Field_t>::point_t Point;
  //   typedef typename std_types<Field_t>::direction_t Direction;
  //   static const int DIM = std_types<Field_t>::n_dim;
        
  //   // fermion
  //   typedef SpinColor<DIM> Fermion;
  //   typedef typename fields::LocalField<Fermion, DIM> FermionField;

  //   // checker board hyper cube size
  //   // c.f. geometry and localfield for more info
  //   static const int n_cb = 0;
    
  //   WilsonBase(const Field_t& G, const FermionField& X, 
  // 			   const double& m ) :
  //     U(G), src(X), mass(m) { };
    
  //   void operator() ( FermionField& dest, const Point& n) {
  //     M(src[n],dest[n],mass);
  //     do_spatial(dest,n);
  //     do_temporal(dest,n,mode_selektor<boundary>());
  //   }
  
  // private:
    
  //   const double& mass;
  //   const Field_t& U;
  //   const FermionField& src;
  //   Mass_t<FermionField> M;
  //   Gamma_t<Field_t,FermionField> gmm;

  //   template <int M> struct mode_selektor { };
  
  //   void do_spatial(FermionField& dest,const Point& n) {
  //     for( Direction mu(1); mu.is_good(); ++mu ) {
  // 	gmm(src,n,mu);
  // 	dest[n] -= (dag(U[n-mu][mu])*gmm.xi1() +
  // 		    U[n ][mu]*gmm.xi2())*.5;
  //     }
  //   }
  //   void do_temporal(FermionField& dest,const Point& n,
  //   		     const mode_selektor<bulk>) {
  //     gmm(src,n,Direction(0));
  //     dest[n] -= (dag(U[n-Direction(0)][Direction(0)])*gmm.xi1() +
  //   		  U[n ][Direction(0)]*gmm.xi2())*.5;
  //   }
  //   void do_temporal(FermionField& dest,const Point& n,
  //   		     const mode_selektor<lower>) {
  //     gmm(src,n,Direction(0));
  //     dest[n] += (dag(U[n-Direction(0)][Direction(0)])*gmm.xi1() -
  //   		  U[n ][Direction(0)]*gmm.xi2())*.5;
  //   }
  //   void do_temporal(FermionField& dest,const Point& n,
  //   		     const mode_selektor<upper>) {
  //     gmm(src,n,Direction(0));
  //     dest[n] -= (dag(U[n-Direction(0)][Direction(0)])*gmm.xi1() -
  //   		  U[n ][Direction(0)]*gmm.xi2())*.5;
  //   }  

  // };


  // template<class Field_t, int boundary> 
  // WilsonTLBase<Field_t,boundary,detail::Mass,detail::Gamma> : WilsonTreeLevelKernel;
  
  // template<class Field_t, int boundary> 
  // typename WilsonBase<Field_t,boundary,detail::Mass,detail::Gamma> WilsonKernel;

  // template<class Field_t, int boundary> 
  // typename WilsonTLBase<Field_t,boundary,detail::Mass5,detail::Gamma5> WilsonTreeLevel5Kernel;
  
  // template<class Field_t, int boundary> 
  // typename WilsonBase<Field_t,boundary,detail::Mass5,detail::Gamma5> Wilson5Kernel;


  //////////////////////////////////////////////////////////////////////
  ///
  ///  TreeLevel ( U = bgf() ) application of the Wilson dirac operator 
  template <class Field_t, int boundary>
  struct WilsonTreeLevelKernel : public base_types<Field_t> {
  public:
    // collect info about the field
    typedef typename std_types<Field_t>::point_t Point;
    typedef typename std_types<Field_t>::direction_t Direction;
    static const int DIM = std_types<Field_t>::n_dim;
        
    // fermion
    typedef SpinColor<DIM> Fermion;
    typedef typename fields::LocalField<Fermion, DIM> FermionField;

    // wilson dirac operator stuff
    typedef typename detail::Mass<FermionField> Mass;
    typedef typename detail::Gamma<Field_t,FermionField> Gamma;
  
    // checker board hyper cube size
    // c.f. geometry and localfield for more info
    static const int n_cb = 0;
    
    WilsonTreeLevelKernel(const Field_t& G, const FermionField& X, 
  			   const double& m ) :
      U(G), src(X), mass(m) { };
    
    void operator() ( FermionField& dest, const Point& n) {
      M(src[n],dest[n],mass);
      do_spatial(dest,n);
      do_temporal(dest,n,mode_selektor<boundary>());
    }
  
  private:
    
    const double& mass;
    const Field_t& U;
    const FermionField& src;
    Mass M;
    Gamma gmm;

    template <int M> struct mode_selektor { };
  
    void do_spatial(FermionField& dest,const Point& n) {
      for( Direction mu(1); mu.is_good(); ++mu ) {
  	gmm(src,n,mu);
  	dest[n] -= (dag(U[n-mu][mu].bgf())*gmm.xi1() +
  		    U[n ][mu].bgf()*gmm.xi2())*.5;
      }
    }
    void do_temporal(FermionField& dest,const Point& n,
    		     const mode_selektor<bulk>) {
      gmm(src,n,Direction(0));
      dest[n] -= (dag(U[n-Direction(0)][Direction(0)].bgf())*gmm.xi1() +
    		  U[n ][Direction(0)].bgf()*gmm.xi2())*.5;
    }
    void do_temporal(FermionField& dest,const Point& n,
    		     const mode_selektor<lower>) {
      gmm(src,n,Direction(0));
      dest[n] += (dag(U[n-Direction(0)][Direction(0)].bgf())*gmm.xi1() -
    		  U[n ][Direction(0)].bgf()*gmm.xi2())*.5;
    }
    void do_temporal(FermionField& dest,const Point& n,
    		     const mode_selektor<upper>) {
      gmm(src,n,Direction(0));
      dest[n] -= (dag(U[n-Direction(0)][Direction(0)].bgf())*gmm.xi1() -
    		  U[n ][Direction(0)].bgf()*gmm.xi2())*.5;
    }  
  };
  
  //////////////////////////////////////////////////////////////////////
  ///
  /// Wilson dirac operator 

  template <class Field_t, int boundary>
  struct WilsonKernel : public base_types<Field_t> {
  public:
    // collect info about the field
    typedef typename std_types<Field_t>::point_t Point;
    typedef typename std_types<Field_t>::direction_t Direction;
    static const int DIM = std_types<Field_t>::n_dim;
        
    // fermion
    typedef SpinColor<DIM> Fermion;
    typedef typename fields::LocalField<Fermion, DIM> FermionField;

    // wilson dirac operator stuff
    typedef typename detail::Mass<FermionField> Mass;
    typedef typename detail::Gamma<Field_t,FermionField> Gamma;
  
    // checker board hyper cube size
    // c.f. geometry and localfield for more info
    static const int n_cb = 0;
    
    WilsonKernel(const Field_t& G, const FermionField& X, 
  			   const double& m ) :
      U(G), src(X), mass(m) { };
    
    void operator() ( FermionField& dest, const Point& n) {
      M(src[n],dest[n],mass);
      do_spatial(dest,n);
      do_temporal(dest,n,mode_selektor<boundary>());
    }
  
  private:
    
    const double& mass;
    const Field_t& U;
    const FermionField& src;
    Mass M;
    Gamma gmm;

    template <int M> struct mode_selektor { };
  
    void do_spatial(FermionField& dest,const Point& n) {
      for( Direction mu(1); mu.is_good(); ++mu ) {
  	gmm(src,n,mu);
  	dest[n] -= (dag(U[n-mu][mu])*gmm.xi1()+U[n ][mu]*gmm.xi2())*.5;
      }
    }
    void do_temporal(FermionField& dest,const Point& n,
    		     const mode_selektor<bulk>) {
      gmm(src,n,Direction(0));
      dest[n] -= (dag(U[n-Direction(0)][Direction(0)])*gmm.xi1() +
    		  U[n ][Direction(0)]*gmm.xi2())*.5;
    }
    void do_temporal(FermionField& dest,const Point& n,
    		     const mode_selektor<lower>) {
      gmm(src,n,Direction(0));
      dest[n] += (dag(U[n-Direction(0)][Direction(0)])*gmm.xi1() -
    		  U[n ][Direction(0)]*gmm.xi2())*.5;
    }
    void do_temporal(FermionField& dest,const Point& n,
    		     const mode_selektor<upper>) {
      gmm(src,n,Direction(0));
      dest[n] -= (dag(U[n-Direction(0)][Direction(0)])*gmm.xi1() -
    		  U[n ][Direction(0)]*gmm.xi2())*.5;
    }  
  };




  //////////////////////////////////////////////////////////////////////
  ///
  ///  TreeLevel ( U = bgf() ) application of the Gamma5 Wilson dirac operator 
  
  template <class Field_t, int boundary>
  struct WilsonTreeLevel5Kernel : public base_types<Field_t> {
  public:
    // collect info about the field
    typedef typename std_types<Field_t>::point_t Point;
    typedef typename std_types<Field_t>::direction_t Direction;
    static const int DIM = std_types<Field_t>::n_dim;
        
    // fermion
    typedef SpinColor<DIM> Fermion;
    typedef typename fields::LocalField<Fermion, DIM> FermionField;

    // wilson dirac operator stuff
    typedef typename detail::Mass5<FermionField> Mass;
    typedef typename detail::Gamma5<Field_t,FermionField> Gamma;
  
    // checker board hyper cube size
    // c.f. geometry and localfield for more info
    static const int n_cb = 0;
    
    WilsonTreeLevel5Kernel(const Field_t& G, const FermionField& X, 
  			   const double& m ) :
      U(G), src(X), mass(m) { };
    
    void operator() ( FermionField& dest, const Point& n) {
      M(src[n],dest[n],mass);
      do_spatial(dest,n);
      do_temporal(dest,n,mode_selektor<boundary>());
    }
  
  private:
    
    const double& mass;
    const Field_t& U;
    const FermionField& src;
    Mass M;
    Gamma gmm;

    template <int M> struct mode_selektor { };
  
    void do_spatial(FermionField& dest,const Point& n) {
      for( Direction mu(1); mu.is_good(); ++mu ) {
  	gmm(src,n,mu);
  	dest[n] -= (dag(U[n-mu][mu].bgf())*gmm.xi1() +
  		    U[n ][mu].bgf()*gmm.xi2())*.5;
      }
    }
    void do_temporal(FermionField& dest,const Point& n,
    		     const mode_selektor<bulk>) {
      gmm(src,n,Direction(0));
      dest[n] -= (dag(U[n-Direction(0)][Direction(0)].bgf())*gmm.xi1() +
    		  U[n ][Direction(0)].bgf()*gmm.xi2())*.5;
    }
    void do_temporal(FermionField& dest,const Point& n,
    		     const mode_selektor<lower>) {
      gmm(src,n,Direction(0));
      dest[n] += (dag(U[n-Direction(0)][Direction(0)].bgf())*gmm.xi1() -
    		  U[n ][Direction(0)].bgf()*gmm.xi2())*.5;
    }
    void do_temporal(FermionField& dest,const Point& n,
    		     const mode_selektor<upper>) {
      gmm(src,n,Direction(0));
      dest[n] -= (dag(U[n-Direction(0)][Direction(0)].bgf())*gmm.xi1() -
    		  U[n ][Direction(0)].bgf()*gmm.xi2())*.5;
    }  
  };
  
  //////////////////////////////////////////////////////////////////////
  ///
  /// Gamma5 Wilson dirac operator 
  template <class Field_t, int boundary>
  struct Wilson5Kernel : public base_types<Field_t> {
  public:
    // collect info about the field
    typedef typename std_types<Field_t>::point_t Point;
    typedef typename std_types<Field_t>::direction_t Direction;
    static const int DIM = std_types<Field_t>::n_dim;
    // fermion
    typedef SpinColor<DIM> Fermion;
    typedef typename fields::LocalField<Fermion, DIM> FermionField;
    // wilson dirac operator stuff
    typedef typename detail::Mass5<FermionField> Mass;
    typedef typename detail::Gamma5<Field_t,FermionField> Gamma;
    // checker board hyper cube size
    static const int n_cb = 0;
    
    Wilson5Kernel(const Field_t& G, const FermionField& X, 
  			   const double& m ) :
      U(G), src(X), mass(m) { };
    
    void operator() ( FermionField& dest, const Point& n) {
      M(src[n],dest[n],mass);
      do_spatial(dest,n);
      do_temporal(dest,n,mode_selektor<boundary>());
    }
  private:
    const double& mass;
    const Field_t& U;
    const FermionField& src;
    Mass M;
    Gamma gmm;

    template <int M> struct mode_selektor { };
  
    void do_spatial(FermionField& dest,const Point& n) {
      for( Direction mu(1); mu.is_good(); ++mu ) {
  	gmm(src,n,mu);
  	dest[n] -= (dag(U[n-mu][mu])*gmm.xi1()+U[n ][mu]*gmm.xi2())*.5;
      }
    }
    void do_temporal(FermionField& dest,const Point& n,
    		     const mode_selektor<bulk>) {
      gmm(src,n,Direction(0));
      dest[n] -= (dag(U[n-Direction(0)][Direction(0)])*gmm.xi1() +
    		  U[n ][Direction(0)]*gmm.xi2())*.5;
    }
    void do_temporal(FermionField& dest,const Point& n,
    		     const mode_selektor<lower>) {
      gmm(src,n,Direction(0));
      dest[n] += (dag(U[n-Direction(0)][Direction(0)])*gmm.xi1() -
    		  U[n ][Direction(0)]*gmm.xi2())*.5;
    }
    void do_temporal(FermionField& dest,const Point& n,
    		     const mode_selektor<upper>) {
      gmm(src,n,Direction(0));
      dest[n] -= (dag(U[n-Direction(0)][Direction(0)])*gmm.xi1() -
    		  U[n ][Direction(0)]*gmm.xi2())*.5;
    }  
  };
  
  
  
  
  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////
  /////
  /////  Perturbative application of the staggered dirac operator
  /////  
  /////
  /////  \author Michele Brambilla <mib.mic@gmail.com>
  /////  \date Fri Nov 02 16:23:46 2012
  //
  //template <class Field_t>
  //struct StaggeredPTKernel {
  //  // collect info about the field
  //  typedef typename std_types<Field_t>::ptGluon_t ptGluon;
  //  typedef typename std_types<Field_t>::ptSU3_t ptSU3;
  //  typedef typename std_types<Field_t>::ptsu3_t ptsu3;
  //  typedef typename std_types<Field_t>::bgf_t BGF;
  //  typedef typename std_types<Field_t>::point_t Point;
  //  typedef typename std_types<Field_t>::direction_t Direction;
  //  static const int ORD = std_types<Field_t>::order;
  //  static const int DIM = std_types<Field_t>::n_dim;
  //
  //  // fermion
  //  typedef SpinColor<1> Fermion;
  //  typedef fields::LocalField<Fermion, DIM> FermionField;
  //  typedef typename array_t<double, DIM>::Type double_array_t;    
  //  
  //  // checker board hyper cube size
  //  // c.f. geometry and localfield for more info
  //  static const int n_cb = 0;
  //
  //
  //
  //  StaggeredPTKernel( Field_t& G, FermionField& X, const double& m ) : U(G), src(X), mass(m) { Uorder = 0; std::cout << "Fix eta definition!!!!!\n";};
  //  
  //  void eta(const Point& n, double_array_t& phase) {
  //
  //  }
  //  
  //  void next_Uorder() { Uorder++; }
  //  void reset_Uorder() { Uorder = 0; }
  //
  //  void operator()( FermionField& dest, Point& n) {
  //    eta(n, phase);
  //    for( Direction nu(0); nu.is_good(); ++nu )
  //	dest[n][0] += .5*phase[(int)nu] * ( U[n][nu][Uorder] * src[n+nu][0] - U[n-nu][nu][Uorder] * src[n-nu][0] );
  //    //      mass
  //    dest[n][0] += mass * src[n][0];
  //
  //  }
  //
  //private:
  //  int Uorder;
  //  double mass;
  //  double_array_t phase;
  //  Field_t& U;
  //  FermionField& src;
  //
  //};


}

#endif




