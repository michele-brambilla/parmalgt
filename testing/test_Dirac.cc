#include "gtest/gtest.h"
#include <Dirac.hpp>
#include <boost/numeric/ublas/io.hpp>

typedef Dirac::DiracAlgebra::matrix_t matrix_t;

struct sim {
  double eps;
  sim(double e = std::numeric_limits<double>::epsilon() * 2) : eps(e)
  { }

  bool operator()(const matrix_t& A, const matrix_t& B){

    bool ret_val = ((A.size1() == B.size1())
                    && (A.size2() == B.size2()));
    if (ret_val)
      for (int i = 0; i < A.size1() && ret_val; ++i)
        for (int j = 0; j < A.size2() && ret_val; ++j)
          ret_val &= (A(i,j) - B(i,j)).abs() < eps;
    return ret_val;
  }
};

TEST(Dirac, AlgebraRelations) {
  using namespace Dirac;
  DiracAlgebra gamma;
  matrix_t two(4,4), zero(4,4);
  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < 4; ++j){
      two(i,j) = 0;
      zero(i,j) = 0;
      if (i == j)
        two(i,j) = 2;
    }
  for (int mu = 0; mu < 6; ++mu)
    for (int nu = 0; nu < 6; ++nu)
      if (mu == nu)
        ASSERT_TRUE(sim()(gamma.a_comm(mu, nu), two));
      else if (mu == 4)
        ASSERT_TRUE(sim()(gamma.a_comm(mu, nu), gamma[nu]*two));
      else if (nu == 4); // same as above, just skip
      else
        ASSERT_TRUE(sim()(gamma.a_comm(mu, nu), zero));
}

TEST(Dirac, Hermitian){
  using namespace Dirac;
  DiracAlgebra gamma;
  matrix_t one(4,4);
  one.clear();
  for (int i = 0; i < 4; ++i)
    one(i,i) = 1;
  for (int mu = 0; mu < 6; ++mu){
    ASSERT_TRUE(sim()(gamma[mu]*dag(gamma[mu]), one));
    ASSERT_TRUE(sim()( dag(gamma[mu]) * gamma[mu], one));
  }
}
