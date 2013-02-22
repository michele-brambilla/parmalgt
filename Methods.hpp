#ifndef METHODS_
#define METHODS_

#include <Kernels.hpp>
#include <Background.h>
#include <LocalField.hpp>
#include <Kernels/generic/RungeKutta.hpp>

namespace meth{

  namespace gf {

    static const int GF_MODE = 1;
    
    namespace detail {
      ////////////////////////////////////////////////////////////
      //
      //  Perform gauge fixing at t = 0. This is the corret method for
      //  an Abelian background field.
      //
      //  \tparam    Fld_t The gluon field type.
      //  \date      Thu Feb 21 18:27:43 2013
      //  \author    Dirk Hesse <dirk.hesse@fis.unipr.it>
      template <class Fld_t>
      void bnd_gauge_fixing(Fld_t& U, const double& alpha, bgf::AbelianBgf){
	typedef kernels::GFMeasKernel<Fld_t> GFMeasKernel;
	typedef kernels::GFApplyKernel<Fld_t> GFApplyKernel;
	int L = U.extent(1);
	GFMeasKernel gfm;
	U.apply_on_timeslice(gfm, 0);
	GFApplyKernel gfa(gfm.val, alpha, L);
	U.apply_on_timeslice(gfa, 0);
      }
      ////////////////////////////////////////////////////////////
      //
      //  Perform gauge fixing at t = 0. This is the corret method for
      //  a trivial background field.
      //
      //  \tparam    Fld_t The gluon field type.
      //  \date      Thu Feb 21 18:28:31 2013
      //  \author    Dirk Hesse <dirk.hesse@fis.unipr.it>
      template <class Fld_t>
      void bnd_gauge_fixing(Fld_t& U, const double& alpha, bgf::ScalarBgf){
	typedef kernels::GFMeasKernel<Fld_t> GFMeasKernel;
	typedef kernels::GFApplyKernelTRBG<Fld_t> GFApplyKernel;
	int L = U.extent(1);
	GFMeasKernel gfm;
	U.apply_on_timeslice(gfm, 0);
	GFApplyKernel gfa(gfm.val, alpha, L);
	U.apply_on_timeslice(gfa, 0);
      }
    }
    ////////////////////////////////////////////////////////////
    //
    //  Perform the complete gauge fixing, at the boundary and in the
    //  bulk. The correct method at the boundary (depending on if we
    //  use a trivial or an Abelian background field is chosen
    //  automatically.
    //
    //  \param     U Gluon field.
    //  \param     alpha Gauge fixing parameter.
    //  \tparam    Fld_t The gluon field type.
    //
    //  \date      Thu Feb 21 18:29:25 2013
    //  \author    Dirk Hesse <dirk.hesse@fis.unipr.it>
    template <class Fld_t>
    void sf_gauge_fixing(Fld_t& U, const double& alpha){
      typedef kernels::GaugeFixingKernel<GF_MODE, Fld_t> GaugeFixingKernel;
      detail::bnd_gauge_fixing(U, alpha, typename kernels::std_types<Fld_t>::bgf_t());
      int T = U.extent(0) - 1;
      GaugeFixingKernel gf(alpha);
      for (int t = 1; t < T; ++t)
	U.apply_on_timeslice(gf, t);
    }
  }

  namespace wflow {

    ////////////////////////////////////////////////////////////
    //
    //  Third-order Runge-Kutta Wilson flow.
    //
    //  \bug       CURRENTLY DOES NOT WORK!!!
    //  \param     U Gluon field.
    //  \param     eps Integration steps size.
    //  \tparam    Fld_t Gluon field type.
    //
    //  \date      Thu Feb 21 18:30:56 2013
    //  \author    Dirk Hesse <dirk.hesse@fis.unipr.it>
    template <class Fld_t>
    void RK3_flow(Fld_t& U, const double& eps){
      typedef typename kernels::flow::WF_RK_1<Fld_t, kernels::StapleSqKernel<Fld_t> > wf1_t;
      typedef typename kernels::flow::WF_RK_2<Fld_t, kernels::StapleSqKernel<Fld_t> > wf2_t;
      typedef typename kernels::flow::WF_RK_3<Fld_t, kernels::StapleSqKernel<Fld_t> > wf3_t;
      typedef typename kernels::std_types<Fld_t>::direction_t Direction;
      Fld_t F(U);
      int T = U.extent(0) - 1;
      std::vector<wf1_t> wf1;
      for (Direction mu; mu.is_good(); ++mu)
        wf1.push_back(wf1_t(mu, eps, F));
      U.apply_on_timeslice(wf1[0], 0);
      for (int t = 1; t < T; ++t)
        for (Direction mu; mu.is_good(); ++mu)
  	U.apply_on_timeslice(wf1[mu], t);
      std::vector<wf2_t> wf2;
      for (Direction mu; mu.is_good(); ++mu)
        wf2.push_back(wf2_t(mu, eps, F));
      U.apply_on_timeslice(wf2[0], 0);
      for (int t = 1; t < T; ++t)
  	for (Direction mu; mu.is_good(); ++mu)
  	  U.apply_on_timeslice(wf2[mu], t);
      std::vector<wf3_t> wf3;
      for (Direction mu; mu.is_good(); ++mu)
        wf3.push_back(wf3_t(mu, eps, F));
      U.apply_on_timeslice(wf3[0], 0);
      for (int t = 1; t < T; ++t)
        for (Direction mu; mu.is_good(); ++mu)
  	U.apply_on_timeslice(wf3[mu], t);
    }
    ////////////////////////////////////////////////////////////
    //
    //  Second-order Runge-Kutta Wilson flow.
    //
    //  \bug       CURRENTLY DOES NOT WORK!!!
    //  \param     U Gluon field.
    //  \param     eps Integration steps size.
    //  \tparam    Fld_t Gluon field type.
    //
    //  \date      Thu Feb 21 18:32:17 2013
    //  \author    Dirk Hesse <dirk.hesse@fis.unipr.it>
    template <class Fld_t>
    void RK2_flow(Fld_t& U, const double& eps){
      typedef typename kernels::flow::WF_RK2_1<Fld_t, kernels::StapleSqKernel<Fld_t> > wf1_t;
      typedef typename kernels::flow::WF_RK2_2<Fld_t, kernels::StapleSqKernel<Fld_t> > wf2_t;
      typedef typename kernels::std_types<Fld_t>::direction_t Direction;
      Fld_t F(U);
      int T = U.extent(0) - 1;
      std::vector<wf1_t> wf1;
      for (Direction mu; mu.is_good(); ++mu)
        wf1.push_back(wf1_t(mu, eps, F));
      U.apply_on_timeslice(wf1[0], 0);
      for (int t = 1; t < T; ++t)
        for (Direction mu; mu.is_good(); ++mu)
  	U.apply_on_timeslice(wf1[mu], t);
      std::vector<wf2_t> wf2;
      for (Direction mu; mu.is_good(); ++mu)
        wf2.push_back(wf2_t(mu, eps, F));
      U.apply_on_timeslice(wf2[0], 0);
      for (int t = 1; t < T; ++t)
        for (Direction mu; mu.is_good(); ++mu)
  	U.apply_on_timeslice(wf2[mu], t);
    }
    ////////////////////////////////////////////////////////////
    //
    //  Wilson flow using an Euler scheme as the integrator. Seems to
    //  work well when performing a simultaneous extrapolation to zero
    //  step size of this method and the gauge update.
    //
    //  \param     U Gluon field.
    //  \param     eps Integration steps size.
    //  \tparam    Fld_t Gluon field type.
    //
    //  \date      Thu Feb 21 18:32:31 2013
    //  \author    Dirk Hesse <dirk.hesse@fis.unipr.it>
    template <class Fld_t>
    void euler_flow(Fld_t& U, const double& eps){

      typedef kernels::StapleSqKernel<Fld_t> StK;
      typedef kernels::TrivialPreProcess<Fld_t> PrK;
      typedef kernels::WilFlowKernel <Fld_t, StK, PrK> WilFlowKernel;
      typedef typename kernels::std_types<Fld_t>::direction_t Direction;
      int T = U.extent(0) - 1;
      StK::weights[0] = 1.;
      std::vector<WilFlowKernel> wf;
      for (Direction mu; mu.is_good(); ++mu)
        wf.push_back(WilFlowKernel(mu, eps));
      // for x_0 = 0 update the temporal direction only
      U.apply_on_timeslice(wf[0], 0);
      // for x_0 != 0 update all directions
      for (int t = 1; t < T; ++t)
	for (Direction mu; mu.is_good(); ++mu)
	  U.apply_on_timeslice(wf[mu], t);
    }
  }

  namespace gu {
    namespace detail{ 
      ////////////////////////////////////////////////////////////
      //
      //  Euler scheme gauge update helper. For the gauge update, we
      //  have to generate an initialze random number generators
      //  residing at each lattice site. This is because we may update
      //  in parallel and hence could have conflicting access to the
      //  same random number generator. However, these must be
      //  initialized once and more importantly the field that holds
      //  them must be generated. To do so, we construct this helper
      //  class and use it as a singleton.
      //  
      //  \tparam   Fld_t The field type, usually LocalField of
      //            something.
      //  \date     Thu Feb 21 18:02:40 2013
      //  \author   Dirk Hesse <dirk.hesse@fis.unipr.it>
      template <class Fld_t>
      struct euler_ {
	typedef kernels::StapleSqKernel<Fld_t> StK;
	typedef kernels::TrivialPreProcess<Fld_t> PrK;
	typedef kernels::GaugeUpdateKernel <Fld_t, StK, PrK> GaugeUpdateKernel;
	typedef typename kernels::std_types<Fld_t>::direction_t Direction;
	int T, L;
	std::vector<GaugeUpdateKernel> gu;
	euler_(const Fld_t& U, const double& eps) : T(U.extent(0) - 1), L(U.extent(1)) {
	  StK::weights[0] = 1.;
	  for (Direction mu; mu.is_good(); ++mu)
	    gu.push_back(GaugeUpdateKernel(mu, eps));
	  GaugeUpdateKernel::rands.resize(L*L*L*(T+1));
	  for (int i = 0; i < L*L*L*(T+1); ++i)
	    GaugeUpdateKernel::rands[i].init(rand());
	}
	////////////////////////////////////////////////////////////
	//
	//  Update a given field. This will take care the whole sweep
	//  over the lattice, respecting the proper degrees of freedom
	//  of the Schrödinger functional
	//
	//  \date      Thu Feb 21 18:18:22 2013
	//  \author    Dirk Hesse <dirk.hesse@fis.unipr.it>
	void operator()(Fld_t& U){
	  int T = U.extent(0) - 1;
	  // for x_0 = 0 update the temporal direction only
	  U.apply_on_timeslice(gu[0], 0);
	  // for x_0 != 0 update all directions
	  for (int t = 1; t < T; ++t)
	    for (Direction mu; mu.is_good(); ++mu)
	      U.apply_on_timeslice(gu[mu], t);
	}
      };
      ////////////////////////////////////////////////////////////
      //
      //  Generate random SU(3) fields. Ther will be one field
      //  associated with each direction. These are needed for a
      //  Runge-Kutta style update of the gauge field. Should be used
      //  as in a singleton to avoid consuming excessive amounts of
      //  memory.
      //
      //  \tparam    Fld_t Gluon field type to be used.
      //  \date      Thu Feb 21 18:19:20 2013
      //  \author    Dirk Hesse <dirk.hesse@fis.unipr.it>
      template <class Fld_t>
      struct rand_gen_ {
	typedef typename kernels::std_types<Fld_t>::SU3_t SU3;
	static const int DIM = kernels::std_types<Fld_t>::n_dim;
	typedef fields::LocalField<SU3, DIM> RandField;
	typedef kernels::RSU3Kernel<RandField> RandKernel;
	typedef typename Fld_t::neighbors_t nt;
	std::vector<RandField> R;
	int L, T;

	rand_gen_(const Fld_t& U) : T(U.extent(0) - 1), 
				    L(U.extent(1)) {
	  RandKernel::rands.resize(L*L*L*(T+1));
	  for (int i = 0; i < L*L*L*(T+1); ++i)
	    RandKernel::rands[i].init(rand());
	  for (int k = 0; k < DIM; ++k)
	    R.push_back(RandField(U.extents(), 1, 0, nt()));
	}
	////////////////////////////////////////////////////////////
	//
	//  Fill the Fields with new random numbers. This should be
	//  called before each update.
	//
	//  \date      Thu Feb 21 18:20:29 2013
	//  \author    Dirk Hesse <dirk.hesse@fis.unipr.it>
	void update() {
	  for (int k = 0; k < DIM; ++k)
	    R[k].apply_everywhere(RandKernel());
	}
	////////////////////////////////////////////////////////////
	//
	//  Access the random SU(3) field associated with a given
	//  direction.
	//
	//  \date      Thu Feb 21 18:22:06 2013
	//  \author    Dirk Hesse <dirk.hesse@fis.unipr.it>
	RandField& operator[](int k){ return R[k]; }
      };
      
    }
    ////////////////////////////////////////////////////////////
    //
    //  Perform a gauge update on the SF d.o.f. using the euler
    //  integratio scheme.
    //
    //  \param     U The gluon field to be updated.
    //  \param     eps Integration step size.
    //
    //  \date      Thu Feb 21 18:22:58 2013
    //  \author    Dirk Hesse <dirk.hesse@fis.unipr.it>
    template <class Fld_t>
    void euler_update(Fld_t& U, const double& eps){
      // singleton euler_ update type
      static detail::euler_<Fld_t> f(U, eps);
      f(U);
    }
    
    ////////////////////////////////////////////////////////////
    //
    //  Perform a gauge update on the SF d.o.f. using a third-order
    //  Runge-Kutta scheme.
    //
    //  \bug       DOES CURRENTLY NOT WORK!
    //
    //  \date      Thu Feb 21 18:26:15 2013
    //  \author    Dirk Hesse <dirk.hesse@fis.unipr.it>
    template <class Fld_t>
    void RK3_update(Fld_t& U, const double& eps){
      typedef typename detail::rand_gen_<Fld_t>::RandField RK_t;
      typedef typename kernels::gauge_update::GU_RK_1<Fld_t, kernels::StapleSqKernel<Fld_t>, RK_t > wf1_t;
      typedef typename kernels::gauge_update::GU_RK_2<Fld_t, kernels::StapleSqKernel<Fld_t>, RK_t > wf2_t;
      typedef typename kernels::gauge_update::GU_RK_3<Fld_t, kernels::StapleSqKernel<Fld_t>, RK_t > wf3_t;
      typedef typename kernels::std_types<Fld_t>::direction_t Direction;
      static detail::rand_gen_<Fld_t> R(U);
      R.update();
      Fld_t F(U);
      int T = U.extent(0) - 1;
      std::vector<wf1_t> wf1;
      for (Direction mu; mu.is_good(); ++mu)
	wf1.push_back(wf1_t(mu, eps, F, R[mu]));
      U.apply_on_timeslice(wf1[0], 0);
      for (int t = 1; t < T; ++t)
	for (Direction mu; mu.is_good(); ++mu)
	  U.apply_on_timeslice(wf1[mu], t);
      std::vector<wf2_t> wf2;
      for (Direction mu; mu.is_good(); ++mu)
	wf2.push_back(wf2_t(mu, eps, F, R[mu]));
      U.apply_on_timeslice(wf2[0], 0);
      for (int t = 1; t < T; ++t)
	for (Direction mu; mu.is_good(); ++mu)
	  U.apply_on_timeslice(wf2[mu], t);
      std::vector<wf3_t> wf3;
      for (Direction mu; mu.is_good(); ++mu)
	wf3.push_back(wf3_t(mu, eps, F, R[mu]));
      U.apply_on_timeslice(wf3[0], 0);
      for (int t = 1; t < T; ++t)
	for (Direction mu; mu.is_good(); ++mu)
	  U.apply_on_timeslice(wf3[mu], t);
    }
    ////////////////////////////////////////////////////////////
    //
    //  Perform a gauge update on the SF d.o.f. using a second-order
    //  Runge-Kutta scheme.
    //
    //  \bug       DOES CURRENTLY NOT WORK!
    //
    //  \date      Thu Feb 21 18:27:23 2013
    //  \author    Dirk Hesse <dirk.hesse@fis.unipr.it>
    template <class Fld_t>
    void RK2_update(Fld_t& U, const double& eps){
      typedef typename detail::rand_gen_<Fld_t>::RandField RK_t;
      typedef typename kernels::gauge_update::GU_RK2_1<Fld_t, kernels::StapleSqKernel<Fld_t>, RK_t > wf1_t;
      typedef typename kernels::gauge_update::GU_RK2_2<Fld_t, kernels::StapleSqKernel<Fld_t>, RK_t > wf2_t;
      typedef typename kernels::std_types<Fld_t>::direction_t Direction;
      static detail::rand_gen_<Fld_t> R(U);
      R.update();
      Fld_t F(U);
      int T = U.extent(0) - 1;
      std::vector<wf1_t> wf1;
      for (Direction mu; mu.is_good(); ++mu)
	wf1.push_back(wf1_t(mu, eps, F, R[mu]));
      U.apply_on_timeslice(wf1[0], 0);
      for (int t = 1; t < T; ++t)
	for (Direction mu; mu.is_good(); ++mu)
	  U.apply_on_timeslice(wf1[mu], t);
      std::vector<wf2_t> wf2;
      for (Direction mu; mu.is_good(); ++mu)
	wf2.push_back(wf2_t(mu, eps, F, R[mu]));
      U.apply_on_timeslice(wf2[0], 0);
      for (int t = 1; t < T; ++t)
	for (Direction mu; mu.is_good(); ++mu)
	  U.apply_on_timeslice(wf2[mu], t);
    }
     ////////////////////////////////////////////////////////////
    //
    //  Perform a gauge update on the SF d.o.f. using a first-order
    //  Runge-Kutta (speak: Euler) scheme.
    //
    //  \bug       NOT TESTED!
    //
    //  \date      Thu Feb 22 12:29:21 2013
    //  \author    Dirk Hesse <dirk.hesse@fis.unipr.it>

    template <class Fld_t>
    void RK1_update(Fld_t& U, const double& eps){
      typedef typename detail::rand_gen_<Fld_t>::RandField RK_t;
      typedef typename kernels::gauge_update::GU_RK1<Fld_t, kernels::StapleSqKernel<Fld_t>, RK_t > wf1_t;
      typedef typename kernels::std_types<Fld_t>::direction_t Direction;
      static detail::rand_gen_<Fld_t> R(U);
      R.update();
      Fld_t F(U);
      int T = U.extent(0) - 1;
      std::vector<wf1_t> wf1;
      for (Direction mu; mu.is_good(); ++mu)
	wf1.push_back(wf1_t(mu, eps, F, R[mu]));
      U.apply_on_timeslice(wf1[0], 0);
      for (int t = 1; t < T; ++t)
	for (Direction mu; mu.is_good(); ++mu)
	  U.apply_on_timeslice(wf1[mu], t);
    }
  }
}

#endif
