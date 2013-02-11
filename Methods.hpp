#ifndef METHODS_
#define METHODS_

#include <Kernels.hpp>
#include <Background.h>

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
      int T = U.extent(0);
      GaugeFixingKernel gf(alpha);
      for (int t = 1; t < T; ++t)
	U.apply_on_timeslice(gf, t);
    }
  }
}

#endif
