#include <Kernels.hpp>
#include <MyMath.h>

#ifndef CLOVER_H_
#define CLOVER_H_

namespace clover {
  namespace detail {
    template <class F>
    struct FieldWrapper {
      typedef typename kernels::std_types<F>::ptSU3_t SU3;
      typedef typename kernels::std_types<F>::direction_t Direction;
      typedef typename kernels::std_types<F>::point_t Point;
      static const int DIM = kernels::std_types<F>::n_dim;
      F *fld;
      FieldWrapper(F& f) : fld(&f) { }
      SU3 operator()(const Point& n, const Direction& mu) const {
        if (mu >= DIM) {
	  return dag((*fld)[n - (-mu)][-mu]);
	}
        else return (*fld)[n][mu];
      }
      SU3 operator()(Point n, const Direction& mu, const Direction& nu) const {
	SU3 tmp = (*this)(n, mu);
	n += mu;
	tmp *= (*this)(n, nu);
	n += nu;
	tmp *= (*this)(n, -mu);
	n -= mu;
	tmp *= (*this)(n, -nu);
	return tmp;
        //return (*this)(n, mu) * (*this)(n + mu, nu) * (*this)(n + nu + mu, -mu) * (*this)(n + nu, -nu);
      }
    };
    template <class F> const int FieldWrapper<F>::DIM;
  } // end namespace detail

  template <class F>
  struct clover {
  };
}

#endif
