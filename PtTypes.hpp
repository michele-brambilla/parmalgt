#ifndef _PT_TYPES_H_
#define _PT_TYPES_H_

#include <Types.h>
#include <MyRand.h>
#include <MyMath.h>
#include <Background.h>
#include <PtAlgo.hpp>

namespace ptt {

  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  ///
  ///  PtMatrix defines a series
  ///
  ///  A = g A^(1) + g^2 A^(2) + ... + g^(ORD + 1) A^(ORD+1)
  ///
  ///  equivalently, we have in matrix notation
  /// 
  ///  A = g A[0] + g^2 A[1] + ... + g^(ORD+1) A[ORD]
  ///
  ///  \author Dirk Hesse <herr.dirk.hesse@gmail.com>
  ///  \date Fri Feb  3 12:39:27 2012

  template <int N> class PtMatrix {
  public:
    typedef PtMatrix self_t;
    typedef typename array_t<SU3, N>::Type su3_array_t;
    typedef typename su3_array_t::value_type SU3_t;
    typedef typename su3_array_t::iterator iterator;
    typedef typename su3_array_t::const_iterator const_iterator;

    //////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////
    ///
    ///  Access operators
    ///
    ///  \author Dirk Hesse <herr.dirk.hesse@gmail.com>
    ///  \date Fri Feb  3 12:41:56 2012

    //const SU3_t& operator[](const int& i) const { return array_[i]; }
    //SU3_t& operator[](const int& i) { return array_[i]; }
    const SU3_t& operator[](const int& i) const { return array_.at(i); }
    SU3_t& operator[](const int& i) { return array_.at(i); }

    //////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////
    ///
    ///  Iterators
    ///
    ///  \author Dirk Hesse <herr.dirk.hesse@gmail.com>
    ///  \date Fri Feb  3 12:42:03 2012

    const_iterator begin() const { return array_.begin(); }
    const_iterator end() const { return array_.end(); }
    iterator begin() { return array_.begin(); }
    iterator end() { return array_.end(); }
  

    //////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////
    ///
    ///  Arithmetics
    ///
    ///  For now, we do not overload operators, but just call the 
    ///  methods e.g add_assign etc.
    ///
    ///  \author Dirk Hesse <herr.dirk.hesse@gmail.com>
    ///  \date Fri Feb  3 12:42:12 2012

    self_t& operator+=(const self_t& other) {
      std::for_each(begin(), end(), pta::incr(other));
      return *this;
    }
  
    self_t operator+(const self_t& other) const {
      self_t result(*this);
      return result += other;
    }

    self_t& operator-=(const self_t& other) {
      std::for_each(begin(), end(), pta::decr(other));
      return *this;
    }
  
    self_t operator-(const self_t& other) const {
      self_t result(*this);
      return result -= other;
    }
  
    template <typename T> self_t& operator*=(const T& other) {
      std::for_each(begin(), end(), pta::mul(other));
      return *this;
    }

    self_t& operator*=(const self_t& other) {
      *this = *this * other;
      return *this;
    }

    template <typename T> self_t& operator/=(const T& other) {
      std::for_each(begin(), end(), pta::div(other));
      return *this;
    }

  private:
    su3_array_t array_;
  };
  
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
///
///  Random PtMatrix
///
///  \author Dirk Hesse <herr.dirk.hesse@gmail.com>
///  \date Fri Feb  3 12:43:06 2012

  template <int N> inline PtMatrix<N> get_random_pt_matrix(){
    PtMatrix<N> result;
    static MyRand r(1235431);
#ifdef HAVE_CXX0X
    for (auto& e : result) { e = SU3rand(r); }
#else
    for (int i = 0; i < N; ++i)
      result[i] = SU3rand(r);
#endif
    return result;
  }

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
///
///  The rest of the arithmetic as templates ...
///
///  \author Dirk Hesse <herr.dirk.hesse@gmail.com>
///  \date Fri Feb  3 12:43:47 2012

  template <int N, typename T> 
  inline PtMatrix<N> operator*(const PtMatrix<N>& A, const T& alpha){
    return PtMatrix<N>(A) *= alpha;
  }

  template <int N>
  inline PtMatrix<N> operator*(const PtMatrix<N>& A, const bgf::AbelianBgf& bg){
#ifdef HAVE_CXX0X
    PtMatrix<N> result(A);
    for (auto& e : result) { e = bg.ApplyFromRight(e); }
#else
    PtMatrix<N> result;
    for (int i = 0; i < N; ++i)
      result[i] = bg.ApplyFromRight(A[i]);
#endif
    return result;
  }

  template <int N>
  inline PtMatrix<N> operator*(const PtMatrix<N>& A, const PtMatrix<N>& B){
    PtMatrix<N> result;
    for (int i = 0; i < N - 1; ++i)
      for (int j = 0; j <= i; ++j)
        result[i + 1] += A[j] * B[i - j];
    return result;
  }
  
  template <int N, typename T> 
  inline PtMatrix<N> operator*(const T& alpha, const PtMatrix<N>& A){
    return PtMatrix<N>(A) *= alpha;
  }

  template <int N>
  inline PtMatrix<N> operator*(const bgf::AbelianBgf& bg, const PtMatrix<N>& A){
#ifdef HAVE_CXX0X
    PtMatrix<N> result(A);
    for (auto& e : result) { e = bg.ApplyFromLeft(e); }
#else
    PtMatrix<N> result;
    for (int i = 0; i < N; ++i)
      result[i] = bg.ApplyFromLeft(A[i]);
#endif
    return result;
  }

}; // namespace ptt

#endif // ndef _PT_TYPES_H_
