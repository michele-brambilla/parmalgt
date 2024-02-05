#ifndef METHODS_
#define METHODS_

#include <background/AbelianBackground.hpp>
#include <background/ScalarBackground.hpp>
#include <kernels/Kernels.hpp>
#include <kernels/Project.hpp>
#include <kernels/RungeKutta.hpp>

namespace meth {

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
void bnd_gauge_fixing(Fld_t &U, const double &alpha, bgf::AbelianBgf) {
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
void bnd_gauge_fixing(Fld_t &U, const double &alpha, bgf::ScalarBgf) {
    typedef kernels::GFMeasKernel<Fld_t> GFMeasKernel;
    typedef kernels::GFApplyKernelTRBG<Fld_t> GFApplyKernel;
    int L = U.extent(1);
    GFMeasKernel gfm;
    U.apply_on_timeslice(gfm, 0);
    GFApplyKernel gfa(gfm.val, alpha, L);
    U.apply_on_timeslice(gfa, 0);
}
} // namespace detail
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
void sf_gauge_fixing(Fld_t &U, const double &alpha) {
    typedef kernels::GaugeFixingKernel<GF_MODE, Fld_t> GaugeFixingKernel;
    detail::bnd_gauge_fixing(U, alpha, typename kernels::std_types<Fld_t>::bgf_t());
    int T = U.extent(0) - 1;
    GaugeFixingKernel gf(alpha);
    for (int t = 1; t < T; ++t)
        U.apply_on_timeslice(gf, t);
}
} // namespace gf

namespace wflow {

namespace detail {
template <class Field_t, class KernelVector_t>
void ApplyOnSfDOF(Field_t &U, KernelVector_t &K) {
    typedef typename kernels::std_types<Field_t>::direction_t Direction;
    int T = U.extent(0) - 1;
    U.apply_on_timeslice(K[0], 0);
    for (int t = 1; t < T; ++t)
        for (Direction mu; mu.is_good(); ++mu)
            U.apply_on_timeslice(K[mu], t);
}
} // namespace detail
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
void RK3_flow(Fld_t &U, const double &eps) {
    using wf1_t = typename kernels::flow::WF_RK_1<Fld_t, kernels::StapleSqKernel<Fld_t>>;
    using wf2_t = typename kernels::flow::WF_RK_2<Fld_t, kernels::StapleSqKernel<Fld_t>>;
    using wf3_t = typename kernels::flow::WF_RK_3<Fld_t, kernels::StapleSqKernel<Fld_t>>;
    using ApplyForceKernel = typename kernels::flow::ApplyForceKernel<Fld_t>;
    using Direction = typename kernels::std_types<Fld_t>::direction_t;
    Fld_t F(U);
    // int T = U.extent(0) - 1;
    std::vector<ApplyForceKernel> afk;
    for (Direction mu; mu.is_good(); ++mu)
        afk.push_back(ApplyForceKernel(mu, eps, F));
    // Step 1
    std::vector<wf1_t> wf1;
    for (Direction mu; mu.is_good(); ++mu)
        wf1.push_back(wf1_t(mu, eps, F));
    detail::ApplyOnSfDOF(U, wf1); // measure
    detail::ApplyOnSfDOF(U, afk); // apply
    // Step 2
    std::vector<wf2_t> wf2;
    for (Direction mu; mu.is_good(); ++mu)
        wf2.push_back(wf2_t(mu, eps, F));
    detail::ApplyOnSfDOF(U, wf2); // measure
    detail::ApplyOnSfDOF(U, afk); // apply
    std::vector<wf3_t> wf3;
    for (Direction mu; mu.is_good(); ++mu)
        wf3.push_back(wf3_t(mu, eps, F));
    detail::ApplyOnSfDOF(U, wf3); // measure
    detail::ApplyOnSfDOF(U, afk); // apply
}
////////////////////////////////////////////////////////////
//
//  Second-order Runge-Kutta Wilson flow.
//
//  \warning   Seems to work, but not thoroughly tested.
//  \param     U Gluon field.
//  \param     eps Integration steps size.
//  \tparam    Fld_t Gluon field type.
//
//  \date      Thu Feb 21 18:32:17 2013
//  \author    Dirk Hesse <dirk.hesse@fis.unipr.it>
template <class Fld_t>
void RK2_flow(Fld_t &U, const double &eps) {
    using wf1_t = typename kernels::flow::WF_RK2_1<Fld_t, kernels::StapleSqKernel<Fld_t>>;
    using wf2_t =  typename kernels::flow::WF_RK2_2<Fld_t, kernels::StapleSqKernel<Fld_t>>;
    using Direction = typename kernels::std_types<Fld_t>::direction_t;
    
    Fld_t F(U);
    Fld_t Util(U);
    int T = U.extent(0) - 1;
    std::vector<wf1_t> wf1;
    for (Direction mu; mu.is_good(); ++mu)
        wf1.push_back(wf1_t(mu, eps, F, Util));
    U.apply_on_timeslice(wf1[0], 0);
    for (int t = 1; t < T; ++t)
        for (Direction mu; mu.is_good(); ++mu)
            U.apply_on_timeslice(wf1[mu], t);
    std::vector<wf2_t> wf2;
    for (Direction mu; mu.is_good(); ++mu)
        wf2.push_back(wf2_t(mu, eps, F, Util));
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
void euler_flow(Fld_t &U, const double &eps) {

    typedef kernels::StapleSqKernel<Fld_t> StK;
    typedef kernels::TrivialPreProcess<Fld_t> PrK;
    typedef kernels::WilFlowKernel<Fld_t, StK, PrK> WilFlowKernel;
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
} // namespace wflow

namespace proj {
template <class Fld_t>
void ProjectSUN(Fld_t &U) {
    using  Kernel_t = typename kernels::proj::ProjectSUN<Fld_t>;
    using Direction = typename kernels::std_types<Fld_t>::direction_t;
    // int T = U.extent(0) - 1;
    std::vector<Kernel_t> proj;
    for (Direction mu; mu.is_good(); ++mu)
        proj.push_back(Kernel_t(mu));
    wflow::detail::ApplyOnSfDOF(U, proj); // measure
}
} // namespace proj

namespace gu {
namespace detail {
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
    using StK = kernels::StapleSqKernel<Fld_t>;
    using PrK = kernels::TrivialPreProcess<Fld_t>;
    using GaugeUpdateKernel = kernels::GaugeUpdateKernel<Fld_t, StK, PrK>;
    using Direction = typename kernels::std_types<Fld_t>::direction_t;
    int L;
    int T;
    
    std::vector<GaugeUpdateKernel> gu;
    euler_(const Fld_t &U, const double &eps) : T(U.extent(0) - 1), L(U.extent(1)) {
        StK::weights[0] = 1.;
        for (Direction mu; mu.is_good(); ++mu)
            gu.push_back(GaugeUpdateKernel(mu, eps));
        GaugeUpdateKernel::rands.resize(L * L * L * (T + 1));
        for (int i = 0; i < L * L * L * (T + 1); ++i)
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
    void operator()(Fld_t &U) {
        // int T = U.extent(0) - 1;
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
//  memory. Can be instructed to use ranlux instead of MyRand.
//
//  \warning   THE RANLUX VERISON IS NOT VERY WELL TESTED
//  \tparam    Fld_t Gluon field type to be used.
//  \date      Thu Feb 21 18:19:20 2013
//  \author    Dirk Hesse <dirk.hesse@fis.unipr.it>
template <class Fld_t>
struct rand_gen_ {
    using SU3 = typename kernels::std_types<Fld_t>::SU3_t;
    static const int DIM = kernels::std_types<Fld_t>::n_dim;
    using  RandField = fields::LocalField<SU3, DIM>;
    using  RandKernel = kernels::RSU3Kernel<RandField>;
    using  nt = typename Fld_t::neighbors_t;
    std::vector<RandField> R;
    int L;
    int T;

    explicit rand_gen_(const Fld_t &U) : L(U.extent(1)), T(U.extent(0) - 1) {
        //*/
        // Using RANLUX
        long vol = (L * L * L * (T + 1));
        std::vector<int> seeds(vol);
        for (long i = 0; i < vol; ++i)
            seeds[i] = rand();
        RandKernel::rands = typename RandKernel::rand_vec_t(seeds.begin(), seeds.end());
        /*/
        // Using MyRand
        RandKernel::rands.resize(L*L*L*(T+1));
        for (int i = 0; i < L*L*L*(T+1); ++i)
          RandKernel::rands[i].init(rand());
        //*/
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
    RandField &operator[](int k) { return R[k]; }
};

} // namespace detail
////////////////////////////////////////////////////////////
//
//  Perform a gauge update on the SF d.o.f. using the euler
//  integratio scheme.
//
//  \param     U The gluon field to be updated.
//  \param     eps Integration step size.
//
//  \warning   YOU MIGHT WANT TO USE THE RK1_update METHOD
//             INSTEAD, WHICH ALLOWS YOU TO USE RANLUX AS RANDOM
//             NUMBER GENERATOR.
//  \date      Thu Feb 21 18:22:58 2013
//  \author    Dirk Hesse <dirk.hesse@fis.unipr.it>
template <class Fld_t>
void euler_update(Fld_t &U, const double &eps) {
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
void RK3_update(Fld_t &U, const double &eps) {
    typedef typename detail::rand_gen_<Fld_t>::RandField RK_t;
    typedef typename kernels::gauge_update::GU_RK_1<Fld_t, kernels::StapleSqKernel<Fld_t>, RK_t> wf1_t;
    typedef typename kernels::gauge_update::GU_RK_2<Fld_t, kernels::StapleSqKernel<Fld_t>, RK_t> wf2_t;
    typedef typename kernels::gauge_update::GU_RK_3<Fld_t, kernels::StapleSqKernel<Fld_t>, RK_t> wf3_t;
    typedef typename kernels::std_types<Fld_t>::direction_t Direction;
    static detail::rand_gen_<Fld_t> R(U);
    R.update();
    Fld_t F(U);
    int T = U.extent(0) - 1;
    std::vector<wf1_t> wf1;
    for (Direction mu; mu.is_good(); ++mu)
        wf1.push_back(wf1_t(int(mu), eps, F, R[int(mu)]));
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

// select the correct kernel for the temporal boundary links
template <class Fld_t, bool b>
struct bnd_kernel { // no ct
    typedef typename detail::rand_gen_<Fld_t>::RandField RK_t;
    typedef typename kernels::StapleSqKernel<Fld_t> Staple_k;
    // kernels to be applied in temporal directions at x_0 = 0, x_0 = T - a
    typedef typename kernels::gauge_update::GU_RK2_1<Fld_t, Staple_k, RK_t> wf1t_t;
    typedef typename kernels::gauge_update::GU_RK2_2<Fld_t, Staple_k, RK_t> wf2t_t;
    // kernels to be applided in spatial directions at x_0 = a, x_0 = T-a
    typedef typename kernels::gauge_update::GU_RK2_1<Fld_t, Staple_k, RK_t> wf1s_early_t;
    typedef typename kernels::gauge_update::GU_RK2_1<Fld_t, Staple_k, RK_t> wf1s_late_t;
    typedef typename kernels::gauge_update::GU_RK2_2<Fld_t, Staple_k, RK_t> wf2s_early_t;
    typedef typename kernels::gauge_update::GU_RK2_2<Fld_t, Staple_k, RK_t> wf2s_late_t;
};
template <class Fld_t>
struct bnd_kernel<Fld_t, true> { // with ct
    typedef typename detail::rand_gen_<Fld_t>::RandField RK_t;
    // plaquette for temporal direction
    typedef typename kernels::StapleSqKernelCtOneTemporal<Fld_t> CtStaple;
    typedef typename kernels::StapleSqKernelCtOneSpatialEarly<Fld_t> CtEarlyStaple;
    typedef typename kernels::StapleSqKernelCtOneSpatialLate<Fld_t> CtLateStaple;
    // kernels to be applied in temporal directions at x_0 = 0, x_0 = T - a
    typedef typename kernels::gauge_update::GU_RK2_1<Fld_t, CtStaple, RK_t> wf1t_t;
    typedef typename kernels::gauge_update::GU_RK2_2<Fld_t, CtStaple, RK_t> wf2t_t;
    // kernels to be applided in spatial directions at x_0 = a, x_0 = T-a
    typedef typename kernels::gauge_update::GU_RK2_1<Fld_t, CtEarlyStaple, RK_t> wf1s_early_t;
    typedef typename kernels::gauge_update::GU_RK2_1<Fld_t, CtLateStaple, RK_t> wf1s_late_t;
    typedef typename kernels::gauge_update::GU_RK2_2<Fld_t, CtEarlyStaple, RK_t> wf2s_early_t;
    typedef typename kernels::gauge_update::GU_RK2_2<Fld_t, CtLateStaple, RK_t> wf2s_late_t;
};

////////////////////////////////////////////////////////////
//
//  Perform a gauge update on the SF d.o.f. using a second-order
//  Runge-Kutta scheme.
//
//  \warning   Seems to work but is not very well tested.
//
//  \date      Thu Feb 21 18:27:23 2013
//  \author    Dirk Hesse <dirk.hesse@fis.unipr.it>
template <class Fld_t, bool do_ct>
void RK2_update(Fld_t &U, const double &eps) {
    using RK_t = typename detail::rand_gen_<Fld_t>::RandField;
    using Staple_k = typename kernels::StapleSqKernel<Fld_t>;
    using wf1_t = typename kernels::gauge_update::GU_RK2_1<Fld_t, Staple_k, RK_t>;
    using wf2_t = typename kernels::gauge_update::GU_RK2_2<Fld_t, Staple_k, RK_t>;
    using Direction = typename kernels::std_types<Fld_t>::direction_t;
    static detail::rand_gen_<Fld_t> R(U);
    R.update();
    Fld_t F(U);
    Fld_t Util(U);
    int T = U.extent(0) - 1;
    ///////////////////////////////
    // FIRST STEP
    std::vector<wf1_t> wf1;
    for (Direction mu; mu.is_good(); ++mu)
        wf1.push_back(wf1_t(mu, eps, F, R[int(mu)], Util));
    // t = 0, mu = 0
    typename bnd_kernel<Fld_t, do_ct>::wf1t_t bndk1(Direction(0), eps, F, R[0], Util);
    U.apply_on_timeslice(bndk1, 0);
    // t = a, mu = 0
    U.apply_on_timeslice(wf1[0], 1);
    // t = a, mu = 1,2,..., D-1
    for (Direction mu(1); mu.is_good(); ++mu) {
        typename bnd_kernel<Fld_t, do_ct>::wf1s_early_t tmp(mu, eps, F, R[mu], Util);
        U.apply_on_timeslice(tmp, 1);
    }
    // t = 2*a, ..., T-2
    for (int t = 2; t < T - 1; ++t)
        for (Direction mu; mu.is_good(); ++mu)
            U.apply_on_timeslice(wf1[mu], t);
    // t = T - 1, mu = 1,2,..., D-1
    for (Direction mu(1); mu.is_good(); ++mu) {
        typename bnd_kernel<Fld_t, do_ct>::wf1s_late_t tmp(mu, eps, F, R[mu], Util);
        U.apply_on_timeslice(tmp, T - 1);
    }
    // t = T - 1, mu = 0
    U.apply_on_timeslice(bndk1, T - 1);
    ///////////////////////////////
    // SECOND STEP
    std::vector<wf2_t> wf2;
    for (Direction mu; mu.is_good(); ++mu)
        wf2.push_back(wf2_t(mu, eps, F, R[mu], Util));
    // t = 0, mu = 0
    typename bnd_kernel<Fld_t, do_ct>::wf2t_t bndk2(Direction(0), eps, F, R[0], Util);
    U.apply_on_timeslice(bndk2, 0);
    // t = a, mu = 0
    U.apply_on_timeslice(wf2[0], 1);
    // t = a, mu = 1,2,..., D-1
    for (Direction mu(1); mu.is_good(); ++mu) {
        typename bnd_kernel<Fld_t, do_ct>::wf2s_early_t tmp(mu, eps, F, R[mu], Util);
        U.apply_on_timeslice(tmp, 1);
    }
    // t = 2*a, ..., T-2
    for (int t = 2; t < T - 1; ++t)
        for (Direction mu; mu.is_good(); ++mu)
            U.apply_on_timeslice(wf2[mu], t);
    // t = T - 1, mu = 1,2,..., D-1
    for (Direction mu(1); mu.is_good(); ++mu) {
        typename bnd_kernel<Fld_t, do_ct>::wf2s_late_t tmp(mu, eps, F, R[mu], Util);
        U.apply_on_timeslice(tmp, T - 1);
    }
    // t = T - 1, mu = 0
    U.apply_on_timeslice(bndk2, T - 1);
}
////////////////////////////////////////////////////////////
//
//  Perform a gauge update on the SF d.o.f. using a first-order
//  Runge-Kutta (speak: Euler) scheme.
//
//  \warning    Not thoroughly tested.
//
//  \date      Thu Feb 22 12:29:21 2013
//  \author    Dirk Hesse <dirk.hesse@fis.unipr.it>

template <class Fld_t>
void RK1_update(Fld_t &U, const double &eps) {
    typedef typename detail::rand_gen_<Fld_t>::RandField RK_t;
    typedef typename kernels::gauge_update::GU_RK1<Fld_t, kernels::StapleSqKernel<Fld_t>, RK_t> wf1_t;
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
} // namespace gu
} // namespace meth

#endif
