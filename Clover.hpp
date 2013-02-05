#include <Kernels.hpp>
#include <MyMath.h>

#ifndef CLOVER_H_
#define CLOVER_H_

namespace clover {
  namespace detail {
    template <class F>
    struct FieldWrapper {
      typedef typename kernels::std_types<F>::ptSU3_t ptSU3;
      typedef typename kernels::std_types<F>::direction_t Direction;
      typedef typename kernels::std_types<F>::point_t Point;
      static const int DIM = kernels::std_types<F>::n_dim;
      F const *fld;
      FieldWrapper(const F& f) : fld(&f) { }
      ptSU3 operator()(const Point& n, const Direction& mu) const {
        if (mu >= DIM) {
	  return dag((*fld)[n - (-mu)][-mu]);
	}
        else return (*fld)[n][mu];
      }
      ptSU3 operator()(Point n, const Direction& mu, const Direction& nu) const {
	ptSU3 tmp = (*this)(n, mu);
	n += mu;
	tmp *= (*this)(n, nu);
	n += nu;
	tmp *= (*this)(n, -mu);
	n -= mu;
	tmp *= (*this)(n, -nu);
	return tmp;
      }
    };
    template <class F> const int FieldWrapper<F>::DIM;
  } // end namespace detail

  template <class Fld_t>
  struct F {
    typedef typename kernels::std_types<Fld_t>::ptSU3_t ptSU3;
    typedef typename kernels::std_types<Fld_t>::point_t Point;
    typedef typename kernels::std_types<Fld_t>::direction_t Direction;
    Direction mu_, nu_;

    F(const Direction& mu, const Direction& nu) : mu_(mu), nu_(nu) { }

    ptSU3 operator()(const Fld_t& UU, const Point& n) const {
      return 1./8 * ( U(n, mu_, nu_) - U(n, nu_, mu_)
		      + U(n, nu_, -mu_) - U(n, mu_, -nu_)
		      + U(n, -mu_, -nu_) - U(n, -nu_, -mu_)
		      + U(n, -nu_, mu_) - U(n, -mu_ ,nu_) );
    }
  };

  template <class Fld_t>
  struct E0s {
    typedef typename kernels::std_types<Fld_t>::ptSU3_t ptSU3;
    typedef typename kernels::std_types<Fld_t>::point_t Point;
    typedef typename kernels::std_types<Fld_t>::direction_t Direction;

    ptSU3 operator()(const Fld_t& UU, const Point& n) const{
      ptSU3 result, tmp;
      for (Direction k(1); k.is_good(); ++k)
	for (Direction l(k + 1); l.is_good(); ++l){
	  F<Fld_t> Fkl(k, l);
	  tmp = Fkl(UU, n);
	  result += tmp*tmp;
	}
      return result;
    }
  };

  template <class Fld_t>
  struct E0m {
    typedef typename kernels::std_types<Fld_t>::ptSU3_t ptSU3;
    typedef typename kernels::std_types<Fld_t>::point_t Point;
    typedef typename kernels::std_types<Fld_t>::direction_t Direction;

    ptSU3 operator()(const Fld_t& UU, const Point& n) const{
      Direction t;
      ptSU3 result, tmp;
      for (Direction k(1); k.is_good(); ++k){
	F<Fld_t> F0k(t, k);
	tmp = F0k(UU, n);
	result += tmp*tmp;
      }
      return result;
    }
  };
}

#endif
