#ifndef POINT_H
#define POINT_H

#include <array>
#include <vector>

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
///
///  D-Dimensional lattice points and directions.
///
///  \author Dirk Hesse <herr.dirk.hesse@gmail.com>
///  \date Mon Mar 26 14:51:20 2012
namespace pt {

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
///
///  Direction class.
///
///  \author Dirk Hesse <herr.dirk.hesse@gmail.com>
///  \date Mon Mar 26 14:51:36 2012
template <int DIM>
class Direction {
  public:
    explicit Direction(const int &m = 0) : mu(m) {}
    Direction &operator++() {
        ++mu;
        return *this;
    }
    Direction &operator--() {
        --mu;
        return *this;
    }
    bool is_good() const { return mu>=0 && mu < DIM; }
    operator int() const { return mu; }
    template <typename A, typename B>
    A deref_fwd(B b) const { return b[DIM + mu]; }
    template <typename A, typename B>
    A deref_bkw(B b) const { return b[mu]; }
    Direction operator-() const { return Direction((mu + DIM) % (2 * DIM)); }
    bool operator<(const int &o) const { return mu < o; }
    bool operator>(const int &o) const { return mu > o; }
    static const Direction t;
    static const Direction x;
    static const Direction y;
    static const Direction z;

  private:
    int mu;
};

template <int DIM>
inline Direction<DIM> operator+(const Direction<DIM> &d, const int &i) {
    Direction<DIM> result(i + d);
    return result;
};

template <int DIM>
class Point {
  public:
    using arr_t = std::array<int, 2 * DIM>;
    using vec_t = std::vector<arr_t>;
    using iter_t = typename vec_t::const_iterator;
    Point(int nn, const iter_t &i) : n(nn), L_begin(i) {}
    Point &operator+=(const Direction<DIM> &mu) {
        if (mu >= DIM) return *this -= Direction<DIM>(mu % DIM);
        n = mu.template deref_fwd<const int &,
                                  const arr_t &>(*(L_begin + n));
        return *this;
    }
    Point &operator-=(const Direction<DIM> &mu) {
        if (mu >= DIM) return *this += Direction<DIM>(mu % DIM);
        n = mu.template deref_bkw<const int &,
                                  const arr_t &>(*(L_begin + n));
        return *this;
    }
    Point &operator++() {
        n++;
        return *this;
    }
    operator int() const { return n; }
    template <typename T>
    typename T::const_reference deref(const T &v) const { return v[n]; }
    template <typename T>
    typename T::reference deref(T &v) const { return v[n]; }
    bool operator==(const Point &other) {
        return n == other.n && L_begin == other.L_begin;
    }
    bool operator!=(const Point &other) {
        return !(*this == other);
    }

  private:
    int n;          // the site
    iter_t L_begin; // positions
};

template <int DIM>
inline Point<DIM> operator+(const Point<DIM> &p, const Direction<DIM> &mu) {
    return Point<DIM>(p) += mu;
}

template <int DIM>
inline Point<DIM> operator-(const Point<DIM> &p, const Direction<DIM> &mu) {
    return Point<DIM>(p) -= mu;
}

} // namespace pt

#endif // ifndef POINT_H
