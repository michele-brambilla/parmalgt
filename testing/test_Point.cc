#include <QCDenvNODEpt.h>
#include <gtest/gtest.h>

using namespace pt;

int size[4] = {5,6,7,8};

struct PointTest : public ::testing::Test {

  PointTest() :  L(size), U(&L), n(U.mk_point(1,2,3,4)) { }
  latt L;
  ptGluon_fld U;
  Point n;
};

TEST_F(PointTest, Chase){
  std::vector<int> y(3);
  for (y[0] = 1; y[0] < size[0]-1; ++y[0])
    for (y[1] = 1; y[1] < size[1]-1; ++y[1])
      for (y[2] = 1; y[2] < size[2]-1; ++y[2])
        for (y[3] = 1; y[3] < size[3]-1; ++y[3])
          for (Direction mu(0); mu.good(); ++mu){
            Point m = U.mk_point(y[0], y[1], y[2], y[3]);
            ++y[mu];
            Point mp = U.mk_point(y[0], y[1], y[2], y[3]);
            for (Direction nu(0); nu.good(); ++nu)
              ASSERT_EQ( &U(m + mu, nu), &U(mp, nu) );
            y[mu] -= 2;
            mp = U.mk_point(y[0], y[1], y[2], y[3]);
            for (Direction nu(0); nu.good(); ++nu)
              ASSERT_EQ( &U(m - mu, nu), &U(mp, nu) );
            ++y[mu];
          }
}

            

TEST_F(PointTest, MoveAround){
  Direction t(0), x(1), y(2), z(3);
  ASSERT_EQ(&U(n + t + x - t - x, t), &U(n, t));
}

TEST_F(PointTest, Wrap){
  for (int i = 0; i < 3; ++i){
    Direction mu(i);
    Point m(n);
    for (int j = 0; j < size[i]; ++j, m += mu);
    ASSERT_TRUE(n == m);
  }
}
