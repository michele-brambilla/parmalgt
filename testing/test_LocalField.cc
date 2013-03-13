#include <gtest/gtest.h>
#include <LocalField.hpp>
#include <cstdlib>
#include <Helper.h>

template <class Field_t>
struct randomize {
  randomize(const int& seed) { std::srand(seed); }
  static const int dim = Field_t::dim;
  void operator()(Field_t& F, const pt::Point<dim>& n) const {
    F[n] = rand();
  }
  static const int n_cb = 0;
};

TEST(Operators, AddAssign){
  typedef fields::LocalField<int, 4> intField;
  intField::extents_t e;
  intField::neighbors_t nn;
  std::fill(e.begin(), e.end(), 10);
  intField A(e,1,0,nn), B(e,1,0,nn);
  A.apply_everywhere(randomize<intField>(123));
  B.apply_everywhere(randomize<intField>(431));
  intField C(A);
  A += B;
  for (intField::const_iterator a = A.begin(), e = A.end(), 
         b = B.begin(), c = C.begin(); a != e; ++a, ++b, ++c)
    ASSERT_EQ(*a, *b + *c);
}

TEST(Operators, SubAssign){
  typedef fields::LocalField<int, 4> intField;
  intField::extents_t e;
  intField::neighbors_t nn;
  std::fill(e.begin(), e.end(), 10);
  intField A(e,1,0,nn), B(e,1,0,nn);
  A.apply_everywhere(randomize<intField>(123));
  B.apply_everywhere(randomize<intField>(431));
  intField C(A);
  A -= B;
  for (intField::const_iterator a = A.begin(), e = A.end(), 
         b = B.begin(), c = C.begin(); a != e; ++a, ++b, ++c)
    ASSERT_EQ(*a, *c - *b);
}

TEST(Operators, InnerProduct){
  typedef fields::LocalField<Cplx, 4> CplxField;
  CplxField::extents_t e;
  CplxField::neighbors_t nn;
  std::fill(e.begin(), e.end(), 10);
  CplxField A(e,1,0,nn), B(e,1,0,nn);
  A.apply_everywhere(randomize<CplxField>(123));
  B.apply_everywhere(randomize<CplxField>(431));
  Cplx manual = 0;
  for (CplxField::const_iterator a = A.begin(), e = A.end(), 
         b = B.begin(); a != e; ++a, ++b)
    manual += (*a)*(*b);
  ASSERT_TRUE(Cmp(manual, A*B)(1e-14));
}

TEST(Operators, ScalarProduct){
  typedef fields::LocalField<Cplx, 4> CplxField;
  CplxField::extents_t e;
  CplxField::neighbors_t nn;
  std::fill(e.begin(), e.end(), 10);
  CplxField A(e,1,0,nn);
  A.apply_everywhere(randomize<CplxField>(123));
  Cplx alpha = 5.234;
  CplxField B(A);
  A *= alpha;
  for (CplxField::const_iterator a = A.begin(), e = A.end(), 
         b = B.begin(); a != e; ++a, ++b)
    ASSERT_TRUE(Cmp(*a, *b * alpha)());
}
