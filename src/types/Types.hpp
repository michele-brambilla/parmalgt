#ifndef _TYPES_H
#define _TYPES_H

#include <array>

#include <complex>
#include <include/SUN.hpp>
#include <include/Vector.hpp>

using SU3 = sun::SU<3>;
using Cplx = sun::SU<3>::data_t;
using CVector = sun::Vec<3>;

using three_vec_t = std::array<Cplx, 3>;

template <class C, int n>
struct array_t {
    using Type = std::array<C, n>;
};

template <int N>
struct Norm : public array_t<double, N>::Type {
    using iterator = typename array_t<double, N>::Type::iterator;
    using const_iterator = typename array_t<double, N>::Type::const_iterator;
    iterator begin() { return rep.begin(); };
    iterator end() { return rep.end(); };
    const_iterator begin() const { return rep.begin(); };
    const_iterator end() const { return rep.end(); };

    explicit Norm(typename array_t<double, N>::Type other) : rep(other){};
    Norm &operator+=(Norm &other) {
        for (int i = 0; i < rep.size(); ++i)
            rep[i] += other[i];
        return *this;
    }
    Norm &operator-=(Norm &other) {
        for (int i = 0; i < rep.size(); ++i)
            rep[i] -= other[i];
        return *this;
    }
    double &operator[](const int i) { return rep[i]; }

  private:
    typename array_t<double, N>::Type rep;
};

template <class C, int n>
std::ostream &operator<<(std::ostream &os, const array_t<C, n> &v) {
    typename array_t<C, n>::iterator it = v.begin();
    for (; it != v.end(); ++it)
        os << *it << "\t";
    return os;
}

template <int NC>
std::ostream &operator<<(std::ostream &os, const sun::SU<NC> &v) {
    for (auto row{0}; row < NC; ++row) {
        for (auto col{0}; col < NC; ++col) {
            os << "(" << v(row, col).real() << "," << v(row, col).imag() << ")\t";
        }
    }
    return os;
}

template <int n>
std::ostream &operator<<(std::ostream &os, const Norm<n> &v) {
    typename Norm<n>::iterator it = v.begin();
    for (; it != v.end(); ++it)
        os << *it << "\t";
    return os;
}

#endif
