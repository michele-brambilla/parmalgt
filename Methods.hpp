#ifndef METHODS_
#define METHODS_

#include <Kernels.hpp>
#include <Background.h>
#include <LocalField.hpp>

namespace meth{

  namespace gf {

    static const int GF_MODE = 1;
    
    namespace detail {
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
  template <class Fld_t>
  void RK3_flow(Fld_t& U, const double& eps){
    typedef typename kernels::WF_RK_1<Fld_t, kernels::StapleSqKernel<Fld_t> > wf1_t;
    typedef typename kernels::WF_RK_2<Fld_t, kernels::StapleSqKernel<Fld_t> > wf2_t;
    typedef typename kernels::WF_RK_3<Fld_t, kernels::StapleSqKernel<Fld_t> > wf3_t;
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
  
  template <class Fld_t>
  void RK2_flow(Fld_t& U, const double& eps){
    typedef typename kernels::WF_RK2_1<Fld_t, kernels::StapleSqKernel<Fld_t> > wf1_t;
    typedef typename kernels::WF_RK2_2<Fld_t, kernels::StapleSqKernel<Fld_t> > wf2_t;
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

  namespace gu {
    template <class Fld_t, class RK_t>
    void RK3_update(Fld_t& U, const double& eps, std::vector<RK_t>& R){
      typedef typename kernels::GU_RK_1<Fld_t, kernels::StapleSqKernel<Fld_t>, RK_t > wf1_t;
      typedef typename kernels::GU_RK_2<Fld_t, kernels::StapleSqKernel<Fld_t>, RK_t > wf2_t;
      typedef typename kernels::GU_RK_3<Fld_t, kernels::StapleSqKernel<Fld_t>, RK_t > wf3_t;
      typedef typename kernels::std_types<Fld_t>::direction_t Direction;
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
    template <class Fld_t, class RK_t>
    void RK2_update(Fld_t& U, const double& eps, std::vector<RK_t>& R){
      typedef typename kernels::GU_RK2_1<Fld_t, kernels::StapleSqKernel<Fld_t>, RK_t > wf1_t;
      typedef typename kernels::GU_RK2_2<Fld_t, kernels::StapleSqKernel<Fld_t>, RK_t > wf2_t;
      typedef typename kernels::std_types<Fld_t>::direction_t Direction;
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
  }
}

#endif
