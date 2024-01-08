#ifndef PROJECT_H_
#define PROJECT_H_

#include "Kernels.hpp"


namespace kernels {
namespace proj {

  template <class Fld_t>
  struct ProjectSUN {
    typedef typename kernels::std_types<Fld_t>::ptSU3_t ptSU3;
    typedef typename kernels::std_types<Fld_t>::point_t Point;
    typedef typename kernels::std_types<Fld_t>::direction_t Direction;
    typedef typename kernels::std_types<Fld_t>::bgf_t BGF;
    static const int ORD = kernels::std_types<Fld_t>::order;
    static const int n_cb = 0;

    ProjectSUN(const Direction& mu) : mu(mu) { }

    void operator()(Fld_t& U, const Point& n) {
      ptSU3 tmp;
      BGF Vinv = U[n][mu].bgf().inverse();
      BGF V = U[n][mu].bgf();
      for (int i = 0; i < ORD; ++i)
	tmp[i] = U[n][mu][i]*Vinv;
      tmp = exp<BGF, ORD>(tmp.reH());
      for (int i = 0; i < ORD; ++i)
	U[n][mu][i] = tmp[i] * V;
    }

    Direction mu;
  };
  
} 
}

#endif
