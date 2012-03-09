#ifndef POINT_H
#define POINT_H

#include <vector>

namespace pt {

  class Direction {
  public:
    explicit Direction(const int& m) : mu(m) { }
    Direction& operator++() { ++mu; return *this; }
    bool good() const { return mu < 4; }
    operator int() const { return mu; }
    template <typename A, typename B>
    A deref_fwd(B b) const { return b[5 + mu]; }
    template <typename A, typename B>
    A deref_bkw(B b) const { return b[mu]; }
    static const Direction t;
    static const Direction x;
    static const Direction y;
    static const Direction z;
  private:
    int mu;
  };

  inline Direction operator+(const Direction& d, const int& i){
    Direction result(i + d);
    return result;
  };

  class Point {
  public:
    typedef std::vector<std::vector<int> >::const_iterator iter_t;
    Point(int nn, const iter_t& i) : n(nn), L_begin(i) { }
    Point& operator+=(const Direction& mu){
      n = mu.deref_fwd<const int&, const std::vector<int>&>(*(L_begin + n));
      return *this;
    }
    Point& operator-=(const Direction& mu){
      n = mu.deref_bkw<const int&, const std::vector<int>&>(*(L_begin + n));
      return *this;
    }
    template <typename T>
    const T& deref(T const * const f) const {
      return f[n];
    }
    template <typename T>
    T& deref(T* const f) const {
      return f[n];
    }
    bool operator==(const Point& other){
      return n == other.n && L_begin == other.L_begin;
    }
  private:
    int n; // the site
    iter_t L_begin; // positions
  };
  
  inline Point operator+(const Point& p, const Direction& mu){
    return Point(p) += mu;
  }
  
  inline Point operator-(const Point& p, const Direction& mu){
    return Point(p) -= mu;
  }

} // namespace pt

#endif //ifndef POINT_H
