#ifndef _BACKGROUND_H_
#define _BACKGROUND_H_

#include "BackgroundInterface.hpp"

#include <ranlxd.hpp>

#include <algorithm>
#include <functional>
#include <iostream>
#include <math.h>
#include <numeric>
#include <vector>

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
/// Background field classes.
///
/// Classes to implement trivial and abelian background fields.

namespace bgf {

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
///
///  Abelian \f$SU(3)\f$ background field.
///
///  \author Dirk Hesse <herr.dirk.hesse@gmail.com>
///  \date Tue Sep 27 11:08:30 2011
class AbelianBgf : public BgfBase {
  public:
    using value_type = three_vec_t::value_type;
    //////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////
    ///
    ///  Vector constructor.
    ///
    ///  Construct from a three_vec_t as the diagonal elements.
    ///
    ///  \param v Vector we want to use for initialzing.
    ///  \author Dirk Hesse <herr.dirk.hesse@gmail.com>
    ///  \date Tue Sep 27 11:07:14 2011
    explicit AbelianBgf(const three_vec_t &v) : v_(v) {}
    AbelianBgf() : v_() { std::fill(v_.begin(), v_.end(), 1); }
    Cplx &operator[](const short &s) {
        return v_[s];
    }
    const Cplx &operator[](const short &s) const {
        return v_[s];
    }
    //////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////
    ///
    ///  Iterators
    ///
    ///  \author Dirk Hesse <herr.dirk.hesse@gmail.com>
    ///  \date Fri Jan 27 13:06:31 2012

    using iterator = three_vec_t::iterator;
    using const_iterator = three_vec_t::const_iterator;
    iterator begin() { return v_.begin(); }
    iterator end() { return v_.end(); }
    const_iterator begin() const { return v_.begin(); }
    const_iterator end() const { return v_.end(); }

    //////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////
    ///
    ///  Left multiplication with an \f$SU(3)\f$ matrix.
    ///
    ///  This returns the product \f$ V U \f$, where \f$ V \f$ is the
    ///  background field's value represented by the instance.
    ///
    ///  \param U The \$f SU3 \$f matrix the instance is to be applied
    ///  to.
    ///  \author Dirk Hesse <herr.dirk.hesse@gmail.com>
    ///  \date Tue Sep 27 11:03:50 2011
    SU3 ApplyFromLeft(const SU3 &U) const override {
        SU3 result;
        for (int i = 0; i < SU3::size; ++i)
            for (int j = 0; j < SU3::size; ++j)
                result(i, j) = v_[i] * U(i, j);
        return result;
    }
    SU3 Add(const SU3 &U) const override {
        SU3 result(U);
        for (int i = 0; i < SU3::size; ++i)
            result(i, i) += v_[i];
        return result;
    }
    CVector ApplyFromLeft(const CVector &v) const override {
        CVector result;
        for (int i = 0; i < CVector::size; ++i)
            result[i] = v_[i] * v[i];
        return result;
    }

    double Norm() const override {
        double result = 0;
#ifdef HAVE_CXX0X
        for (const auto &e : v_) {
            result += abs(e);
        }
#else
        for (int i = 0; i < 3; ++i)
            result += abs(v_[i]);
#endif
        return result;
    }

    //////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////
    ///
    ///  Right multiplication with an \f$SU(3)\f$ matrix.
    ///
    ///  This returns the product \f$ U V \f$, where \f$ V \f$ is the
    ///  background field's value represented by the instance.
    ///
    ///  \param U The \$f SU3 \$f matrix the instance is to be applied
    ///  to.
    ///  \author Dirk Hesse <herr.dirk.hesse@gmail.com>
    ///  \date Tue Sep 27 11:03:50 2011
    SU3 ApplyFromRight(const SU3 &U) const override {
        SU3 result;
        for (int i = 0; i < SU3::size; ++i)
            for (int j = 0; j < SU3::size; ++j)
                result(i, j) = v_[j] * U(i, j);
        return result;
    }
    CVector ApplyFromRight(const CVector &v) const override {
        CVector result;
        for (int i = 0; i < CVector::size; ++i)
            result[i] = v_[i] * v[i];
        return result;
    }
    bool operator==(const AbelianBgf &other) const {
        return v_ == other.v_;
    }
    //////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////
    ///
    ///  Generic *= operator template.
    ///
    ///  \author Dirk Hesse <herr.dirk.hesse@gmail.com>
    ///  \date Mon Sep 26 18:28:23 2011
    template <class C>
    AbelianBgf &operator*=(const C &alpha) {
        for(auto &element : v_) element *= alpha;
        return *this;
    }
    //////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////
    ///
    ///  Generic /= operator template.
    ///
    ///  \author Dirk Hesse <herr.dirk.hesse@gmail.com>
    ///  \date Mon Sep 26 18:28:23 2011
    template <class C>
    AbelianBgf &operator/=(const C &alpha) {
        *this *= 1. / alpha;
        return *this;
    }
    //////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////
    ///
    ///  Generic multiplication.
    ///
    ///  \author Dirk Hesse <herr.dirk.hesse@gmail.com>
    ///  \date Mon Sep 26 18:28:36 2011
    template <class C>
    AbelianBgf operator*(const C &alpha) const {
        AbelianBgf result(*this);
        return result *= alpha;
    }

    // MB: render using free functions
    // //////////////////////////////////////////////////////////////////////
    // //////////////////////////////////////////////////////////////////////
    // ///
    // ///  SpinColor multiplication. Ok, it's awful.. let's see if it works!
    // ///
    // ///  \author Michele Brambilla <mib.mic@gmail.com>
    // ///  \date Fri Feb 01 20:23:23 2013
    // template<int DIM>
    //   SpinColor<DIM> operator* (const SpinColor<DIM>& alpha ) const {
    //   SpinColor<DIM> result(alpha);
    //   for( int mu = 0; mu < DIM; ++mu)
    // 	for( int col = 0; col < 3; ++col)
    // 	  result[mu][col] *= v_[col];
    //   return result;
    // }

    //////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////
    ///
    ///  Addition and subtraction of a sclar.
    ///
    ///  \author Dirk Hesse <herr.dirk.hesse@gmail.com>
    ///  \date Wed Jan 11 18:03:06 2012
    template <class C>
    AbelianBgf &operator+=(const C &alpha) {
        std::transform(v_.begin(), v_.end(), v_.begin(),
                    //    std::bind1st(std::plus<C>(), alpha)
                        [alpha](three_vec_t::value_type &value) { value += alpha; }
                        );

        return *this;
    }
    AbelianBgf &operator+=(const AbelianBgf &alpha) {
        std::transform(v_.begin(), v_.end(), alpha.v_.begin(),
                       v_.begin(), std::plus<Cplx>());
        return *this;
    }
    AbelianBgf operator-() const {
        AbelianBgf result;
        for (int i = 0; i < 3; ++i)
            result[i] = -v_[i];
        return result;
    }
    template <class C>
    AbelianBgf &operator-=(const C &alpha) {
        return *this += (-alpha);
    }
    template <class C>
    AbelianBgf operator+(const C &alpha) const {
        AbelianBgf result(*this);
        return result += alpha;
    }
    template <class C>
    AbelianBgf operator-(const C &alpha) const {
        AbelianBgf result(*this);
        return result -= alpha;
    }
    AbelianBgf inverse() const {
        AbelianBgf result;
        for (int i = 0; i < 3; ++i)
            result[i] = 1. / v_[i];
        return result;
    }
    void set_to_one() override { std::fill(v_.begin(), v_.end(), 1); }
    void set_to_zero() override { std::fill(v_.begin(), v_.end(), 0); }
    /// Trace
    Cplx Tr() const override {
        return std::accumulate(v_.begin(), v_.end(), Cplx(0));
    }

    /// Make traceless
    void Trless() override {
        Cplx Tro3 = Tr() / 3.;
#ifdef HAVE_CXX0X
        for (auto &e : v_)
            e -= Tro3;
#else
        for (int i = 0; i < 3; ++i)
            v_[i] -= Tro3;
#endif
    }

    //////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////
    ///
    ///  reH
    ///
    ///  This takes the traceless part of 0.5*[V - V^\dagger]
    ///
    ///  \author Dirk Hesse <herr.dirk.hesse@gmail.com>
    ///  \date Fri Feb  3 16:23:25 2012

    void reH() override {
#ifdef HAVE_CXX0X
        for (auto &e : v_)
            e.real() = 0;
#else
        for (int i = 0; i < 3; ++i)
            v_[i].real() = 0;
#endif
        Trless();
    }

    AbelianBgf dag() const {
        AbelianBgf result(*this);
        for (int i = 0; i < 3; ++i)
            result[i] = conj(result[i]);
        return result;
    }

    //////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////
    ///
    ///  Number of doubles needed to store the object.
    ///
    ///  \author Dirk Hesse <herr.dirk.hesse@gmail.com>
    ///  \date Sun Mar 25 14:50:37 2012

    static const int storage_size = 6;

    //////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////
    ///
    ///  Buffer to a vector of doubles.
    ///
    ///  \author Dirk Hesse <herr.dirk.hesse@gmail.com>
    ///  \date Sun Mar 25 14:50:53 2012

    std::vector<double>::iterator &
    buffer(std::vector<double>::iterator &i) const {
        for (const_iterator n = begin(); n != end(); ++n) {
            *i = n->real();
            ++i;
            *i = n->imag();
            ++i;
        }
        return i;
    }

    //////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////
    ///
    ///  Write to a file
    ///
    ///  \author Dirk Hesse <herr.dirk.hesse@gmail.com>
    ///  \date Fri May 25 16:02:33 2012

    template <class Writer_t>
    void write(Writer_t &o) const {
        for (auto n = begin(); n != end(); ++n)
            o.write(*n);
    }
    template <class Reader_t>
    void read(Reader_t &o) {
        for (auto n = begin(); n != end(); ++n)
            o.read(*n);
    }

    //////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////
    ///
    ///  Read from a buffer
    ///
    ///  \author Dirk Hesse <herr.dirk.hesse@gmail.com>
    ///  \date Mon Mar 26 16:43:22 2012

    std::vector<double>::const_iterator &
    unbuffer(std::vector<double>::const_iterator &i) {
        for (auto n = begin(); n != end(); ++n) {
            n->real() = *i;
            ++i;
            n->imag() = *i;
            ++i;
        }
        return i;
    }

  private:
    //////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////
    ///
    ///  Diagonal elements of V.
    ///
    ///  \author Dirk Hesse <herr.dirk.hesse@gmail.com>
    ///  \date Tue Sep 27 10:50:29 2011
    three_vec_t v_;
};

inline AbelianBgf dag(const AbelianBgf &b) {
    return b.dag();
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
///
///  Specialization of *= and * operator for AbelianBgf type.
///
///  The Arithmetics are quite straight forward, altough the
///  AbelianBgf*AbelianBgf multiplication is a special case.
///
///  \author Dirk Hesse <herr.dirk.hesse@gmail.com>
///  \date Mon Sep 26 18:25:58 2011

template <>
inline AbelianBgf &AbelianBgf::operator*=(const AbelianBgf &other) {
    std::transform(v_.begin(), v_.end(), other.v_.begin(),
                   v_.begin(), std::multiplies<Cplx>());
    return *this;
}

inline std::ostream &operator<<(std::ostream &os, const AbelianBgf &b) {
    os << "{ ";
    for (auto i = b.begin();
         i != b.end(); ++i)
        os << "(" << i->real() << ", " << i->imag() << ") ";
    os << "}";
    return os;
}

class AbelianBgfFactory {
  public:
    //////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////
    ///
    ///  Initialization of Vt_.
    ///
    ///  Here, we construct \f$ V(t) \f$ a la Peter Weisz, viz.
    ///  \f[ V(t) = \exp i ( {\cal E} x_0 - iC) \f],
    ///  where
    ///  \f[ {\cal E} = -i (C' - C) / T, \f]
    ///  and $C',C$ are the gluon boundary fields.
    ///
    ///  \param T Temporal lattice extend.
    ///  \param L Spatial lattice extend.
    ///  \param eta Background field parameter. C.f. P. Weisz. Usually
    ///  0 is used.
    ///  \param nu Background field parameter. C.f. P. Weisz. Usually
    ///  0 is used.
    ///
    ///  \author Dirk Hesse <herr.dirk.hesse@gmail.com>
    ///  \date Tue Sep 27 10:46:16 2011
    AbelianBgfFactory(const int &T_in, const int &L, const int &s = 0)
        : Vt_(T_in + 1), T(T_in) {
        if (T < 0)
            throw std::exception();
        // if (!s) {
        //   double pi = std::atan(1.)*4.;
        //   double gamma = 1./L/T * (eta + pi/3.);
        //   double eps[3] = {-2.*gamma, gamma, gamma};
        //   double iC[3] = {-( eta - pi/3) / L,
        //                   -eta * (-0.5 + nu)  / L,
        //                   -( eta * (0.5 + nu) + pi/3.) / L};
        //   for (int t = 0; t <= T; ++t)
        //     for (int k = 0; k < 3; ++k)
        //       Vt_[t][k] = exp(Cplx(0,eps[k] * t - iC[k]));
        // }
        //  tree level c_t
        if (T + s != L)
            std::cout
                << "WARNING!\nYou gave strange values for T,L, and s!\n"
                << "T = " << T << "\nL = " << L << "\ns = " << s
                << "\nYou should have T + s = L!\n"
                << "You're entering a world of pain!\n"
                << "This is my last warning.\n";

        double ct = 2. / (2 + s);
        // eta and nu remain zero since we only plugged in f for these values!
        const double eta = 0;
        const double nu = 0;
        double pi = std::atan(1.) * 4.;
        SU3 C;
        SU3 Cp;
        SU3 B;
        C(0, 0) = eta - pi / 3;
        C(1, 1) = eta * (nu - 0.5);
        C(2, 2) = -eta * (nu + 0.5) + pi / 3;
        Cp(0, 0) = -eta - pi;
        Cp(1, 1) = eta * (nu + 0.5) + pi / 3;
        Cp(2, 2) = -eta * (nu - 0.5) + 2. * pi / 3;
        C *= Cplx(0, 1. / L);
        Cp *= Cplx(0, 1. / L);
        // for (int t = 0; t <= T; ++t)
        //   for (int k = 0; k < 3; ++k)
        //     Vt_[t][k] = exp((t * Cp(k,k) + (L - t)*C(k,k))/L);
        //  Values for the diagonals of V at t = 0, T
        for (int k = 0; k < 3; ++k) {
            Vt_[0][k] = exp(C(k, k)) * ct;
            Vt_[T][k] = exp(Cp(k, k)) * ct;
        }
        auto i = Cplx(0, 1);
        for (int t = 1; t < T; ++t) {
            Vt_[t][0] = exp(-2. * f[s + 1][L / 2 - 2] * (t - 0.5 * T) * i + 0.5 * (C(0, 0) + Cp(0, 0)));
            for (int k = 1; k < 3; ++k)
                Vt_[t][k] = exp(f[s + 1][L / 2 - 2] * (t - 0.5 * T) * i + 0.5 * (C(k, k) + Cp(k, k)));
        }
    };
    // spatial components, V_k
    const three_vec_t &get(int t) { return Vt_.at(t); }
    // temporal, V_0, set V_0(t = T) to 0, such that one can exploit
    // the periodicity of the lattice for e.g. LW improved gauge
    // actions.
    three_vec_t get_temporal(int t) const {
        static three_vec_t one = {1, 1, 1};
        static three_vec_t zero = {0, 0, 0};
        if (t >= T)
            return zero;
        return one;
    }

  private:
    std::vector<three_vec_t> Vt_;
    static const std::array<std::array<double, 31>,3> f;
    // static const double f[3][31];
    int T;
};

/// Returns V_\mu(t)
inline AbelianBgf get_abelian_bgf(const int &t,
                                  const int &mu,
                                  const int &T = -1,
                                  const int &L = 0,
                                  const int &s = 0) {
    static AbelianBgfFactory factory(T, L, s);

    if (!mu)
        return AbelianBgf(factory.get_temporal(t));
    else
        return AbelianBgf(factory.get(t));
};

template <class B>
inline B unit() { return B(); }

template <>
inline AbelianBgf unit<AbelianBgf>() {
    // three_vec_t alpha_v = {1, 1, 1};
    return AbelianBgf{three_vec_t{1, 1, 1}};
};

template <class B>
inline B zero() { return B(); }

template <>
inline AbelianBgf zero<AbelianBgf>() {
    three_vec_t alpha_v = {0, 0, 0};
    return AbelianBgf(alpha_v);
};

inline AbelianBgf random() {
    static ranlxd::Rand r(134);
    std::array<double, 2 * 3> storage;
    r.ranlxd(storage.begin(), storage.end());
    three_vec_t alpha_v{Cplx{storage[0], storage[1]},
                        Cplx{storage[2], storage[3]},
                        Cplx{storage[4], storage[5]}};

    return AbelianBgf(alpha_v);
};
} // namespace bgf

#endif
