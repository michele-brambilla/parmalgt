#ifndef _SCALAR_BACKGROUND_HPP_
#define _SCALAR_BACKGROUND_HPP_

#include "BackgroundInterface.hpp"

#include <include/SUN.hpp>
#include <types/Types.hpp>


namespace bgf {

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
///
///  Scalar \f$SU(3)\f$ background field.
///
///  \author Dirk Hesse <herr.dirk.hesse@gmail.com>
///  \date Tue Sep 27 11:08:30 2011
class ScalarBgf : public BgfBase {
  public:
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
    explicit ScalarBgf(const Cplx &a) : a_(a) {}
    ScalarBgf() : a_(1, 0) {}

    virtual SU3 ApplyFromLeft(const SU3 &U) const {
        SU3 result;
        for (int i = 0; i < SU3::rep_size; ++i)
            result[i] = a_ * U[i];
        return result;
    }
    virtual SU3 Add(const SU3 &U) const {
        SU3 result(U);
        for (int i = 0; i < SU3::size; ++i)
            result(i, i) += a_;
        return result;
    }
    virtual CVector ApplyFromLeft(const CVector &v) const {
        CVector result;
        for (int i = 0; i < CVector::size; ++i)
            result[i] = a_ * v[i];
        return result;
    }

    double Norm() const {
        return abs(a_);
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
    virtual SU3 ApplyFromRight(const SU3 &U) const {
        SU3 result;
        for (int i = 0; i < SU3::rep_size; ++i)
            result[i] = a_ * U[i];
        return result;
    }
    virtual CVector ApplyFromRight(const CVector &v) const {
        CVector result;
        for (int i = 0; i < CVector::rep_size; ++i)
            result[i] = a_ * v[i];
        return result;
    }
    bool operator==(const ScalarBgf &other) const {
        return a_ == other.a_;
    }
    //////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////
    ///
    ///  Generic *= operator template.
    ///
    ///  \author Dirk Hesse <herr.dirk.hesse@gmail.com>
    ///  \date Mon Sep 26 18:28:23 2011
    template <class C>
    ScalarBgf &operator*=(const C &alpha) {
        a_ *= alpha;
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
    ScalarBgf &operator/=(const C &alpha) {
        a_ /= alpha;
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
    ScalarBgf operator*(const C &alpha) const {
        ScalarBgf result(*this);
        return result *= alpha;
    }

// MB
//     //////////////////////////////////////////////////////////////////////
//     //////////////////////////////////////////////////////////////////////
//     ///
//     ///  Multiplication with CVector and SpinColor
//     ///
//     ///  \author Michele
//     ///  \date Wed Oct 17
//     // inline CVector operator* (const CVector& v ) const {
//     //   CVector result;
//     //   return result = v * a_;
//     // }
//     template <int DIM>
//     inline SpinColor<DIM> operator*(const SpinColor<DIM> &v) const {
//         SpinColor<DIM> result;
// #ifdef HAVE_CXX0X
//         for (auto it = v.begin(), r = result.begin();
//              it != v.end(); ++it, ++r)
//             (*r) = a_ * (*it);
// #else
//         typename SpinColor<DIM>::const_iterator i = v.begin(), e = v.end();
//         typename SpinColor<DIM>::iterator r = result.begin();
//         for (; i != e; ++i, ++r)
//             (*r) = a_ * (*i);
// #endif
//         return result;
//     }

    //////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////
    ///
    ///  Multiplication with self type.
    ///
    ///  \author Dirk Hesse <herr.dirk.hesse@gmail.com>
    ///  \date Tue Sep  4 15:18:20 2012

    inline ScalarBgf operator*(const ScalarBgf &alpha) const {
        ScalarBgf result(*this);
        return result *= alpha.a_;
    }
    //////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////
    ///
    ///  Addition and subtraction of a sclar.
    ///
    ///  \author Dirk Hesse <herr.dirk.hesse@gmail.com>
    ///  \date Wed Jan 11 18:03:06 2012

    template <class C>
    ScalarBgf &operator+=(const C &alpha) {
        a_ += alpha;
        return *this;
    }
    ScalarBgf &operator+=(const ScalarBgf &alpha) {
        a_ += alpha.a_;
        return *this;
    }
    ScalarBgf operator-() const {
        ScalarBgf result;
        result.a_ = -a_;
        return result;
    }
    template <class C>
    ScalarBgf &operator-=(const C &alpha) {
        a_ -= alpha;
        return *this;
    }
    ScalarBgf &operator-=(const ScalarBgf &alpha) {
        a_ -= alpha.a_;
        return *this;
    }

    template <class C>
    ScalarBgf operator+(const C &alpha) const {
        ScalarBgf result(*this);
        return result += alpha;
    }
    template <class C>
    ScalarBgf operator-(const C &alpha) const {
        ScalarBgf result(*this);
        return result -= alpha;
    }
    ScalarBgf inverse() const {
        ScalarBgf result;
        result.a_ = 1. / a_;
        return result;
    }
    void set_to_one() { a_ = Cplx(1, 0); };
    void set_to_zero() { a_ = Cplx(0, 0); };
    /// Trace
    Cplx Tr() const {
        return a_ * 3.;
    }

    /// Make traceless
    virtual void Trless() {
        set_to_zero();
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

    virtual void reH() {
        Trless();
    }

    ScalarBgf dag() const {
        ScalarBgf result(conj(a_));
        return result;
    }

    //////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////
    ///
    ///  Number of doubles needed to store the object.
    ///
    ///  \author Dirk Hesse <herr.dirk.hesse@gmail.com>
    ///  \date Sun Mar 25 14:50:37 2012

    static const int storage_size = 2;

    //////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////
    ///
    ///  Buffer to a vector of doubles.
    ///
    ///  \author Dirk Hesse <herr.dirk.hesse@gmail.com>
    ///  \date Sun Mar 25 14:50:53 2012

    std::vector<double>::iterator &
    buffer(std::vector<double>::iterator &i) const {
        *i = a_.real();
        ++i;
        *i = a_.imag();
        ++i;
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
        o.write(a_);
    }
    template <class Reader_t>
    void read(Reader_t &o) {
        o.read(a_);
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
        a_.real() = *i;
        ++i;
        a_.imag() = *i;
        ++i;
        return i;
    }

    const Cplx &val() const { return a_; }

  private:
    //////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////
    ///
    ///  Diagonal elements of V.
    ///
    ///  \author Dirk Hesse <herr.dirk.hesse@gmail.com>
    ///  \date Tue Sep 27 10:50:29 2011
    Cplx a_;
};

inline ScalarBgf dag(const ScalarBgf &b) {
    return b.dag();
}

inline std::ostream &operator<<(std::ostream &os, const ScalarBgf &b) {
    os << "{ ";
    os << "(" << b.val().real() << ", " << b.val().imag() << ") ";
    os << "}";
    return os;
}

template <>
inline ScalarBgf zero<ScalarBgf>() {
    return ScalarBgf(0);
};

} // namespace bgf

#endif // _SCALAR_BACKGROUND_HPP_