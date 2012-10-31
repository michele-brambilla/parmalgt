#include "Background.h"
#include "gtest/gtest.h"
#include "Point.hpp"
#include "Geometry.hpp"

const int DIM = 4;
const int SIZE = 8;


TEST(Geometry, Neighbors){

  geometry::Geometry<DIM>::extents_t e;
  std::fill(e.begin(), e.end(), SIZE);
  geometry::Geometry<DIM> g(e);
  geometry::Geometry<DIM>::raw_pt_t n;
  std::fill(n.begin(), n.end(), 0);
  for (; n[0] < SIZE; ++n[0])
    for (; n[1] < SIZE; ++n[1])
      for (; n[2] < SIZE; ++n[2])
        for (; n[3] < SIZE; ++n[3])
          for (pt::Direction<DIM> mu; mu.is_good(); ++mu){
            geometry::Geometry<DIM>::raw_pt_t x(n);
            x[mu]++;
            ASSERT_EQ(g.mk_point(x), g.mk_point(n) + mu);
            ASSERT_NE(g.mk_point(x), g.mk_point(n));
            x[mu] -= 2;
            ASSERT_EQ(g.mk_point(x), g.mk_point(n) - mu);
            ASSERT_NE(g.mk_point(x), g.mk_point(n));
        }
};

TEST(Geometry, ManualVsDirectionAccess){
  geometry::Geometry<DIM>::extents_t e;
  std::fill(e.begin(), e.end(), SIZE);
  geometry::Geometry<DIM> g(e);
  geometry::Geometry<DIM>::raw_pt_t n;
  std::fill(n.begin(), n.end(), 0);
  pt::Point<DIM> x = g.mk_point(n);
  for (; n[0] < SIZE; ++n[0], x += pt::Direction<DIM>(0))
    for (; n[1] < SIZE; ++n[1], x += pt::Direction<DIM>(1))
      for (; n[2] < SIZE; ++n[2], x += pt::Direction<DIM>(2))
        for (; n[3] < SIZE; ++n[3], x += pt::Direction<DIM>(3))
            ASSERT_EQ(x, g.mk_point(n));

};


TEST(Geometry, BinFunction){
  geometry::Geometry<DIM>::extents_t e;
  std::fill(e.begin(), e.end(), SIZE);
  e[1] += 10;
  e[2] += 5;
  e[3] += 12;
  geometry::Geometry<DIM> g(e);
  geometry::Geometry<DIM>::raw_pt_t n;
  std::fill(n.begin(), n.end(), 0);
  pt::Point<DIM> x = g.mk_point(n);
  for (; n[0] < SIZE; ++n[0], x += pt::Direction<DIM>(0))
    for (; n[1] < SIZE; ++n[1], x += pt::Direction<DIM>(1))
      for (; n[2] < SIZE; ++n[2], x += pt::Direction<DIM>(2))
        for (; n[3] < SIZE; ++n[3], x += pt::Direction<DIM>(3)){
          for (pt::Direction<DIM> mu; mu.is_good(); ++mu){
            ASSERT_NE(g.bin(x,1), g.bin(x + mu, 1));
            ASSERT_EQ(g.bin(x,1), g.bin(x + mu + mu, 1));
            for (pt::Direction<DIM> nu; nu.is_good(); ++nu){
              ASSERT_NE(g.bin(x,2), g.bin(x + mu + nu, 2));
              for (pt::Direction<DIM> rho; rho.is_good(); ++rho)
                ASSERT_NE(g.bin(x,2), g.bin(x + mu + nu + rho, 2));
            }
            ASSERT_EQ(g.bin(x,2), g.bin(x + mu + mu + mu + mu, 2));
          }
        }
}
TEST(Geometry, CheckerBoard){
  geometry::Geometry<DIM>::extents_t e;
  std::fill(e.begin(), e.end(), SIZE);
  e[1] += 10;
  e[2] += 5;
  e[3] += 12;
  geometry::Geometry<DIM> g(e);
  typedef geometry::CheckerBoard<DIM> cb_t;
  cb_t cb[3] = {cb_t(g, 1), cb_t(g, 2), cb_t(g, 3)};
  typedef geometry::CheckerBoard<DIM>::slice slice;
  typedef geometry::CheckerBoard<DIM>::bin bin;
  geometry::Geometry<DIM>::raw_pt_t n;
  std::fill(n.begin(), n.end(), 0);
  for (n[0] = 0; n[0] < e[0]; ++n[0])
    for (n[1] = 0; n[1] < e[1]; ++n[1])
      for (n[2] = 0; n[2] < e[2]; ++n[2])
        for (n[3] = 0; n[3] < e[3]; ++n[3]){
          pt::Point<DIM> x = g.mk_point(n);
          for (int k = 0; k < 3; ++k){
            bool found = false;
            // note: use nca instead of operator[] below for
            // non-constant access which we will use to erase the found
            // elements to assure we don't have duplicates or somehting
            // the like (which would be very strange!)
            for (slice::iterator i = cb[k].nca(n[0]).begin();
                 i != cb[k].nca(n[0]).end() && !found; ++i){
              bin::iterator f = std::find(i->begin(), i->end(), x);
              if (f != i->end()){
                found = true;
                i->erase(f);
              }
            }
            ASSERT_EQ(found, true);
          }
        }
  for (int k = 0; k < 3; ++k)
    for (int t = 0; t < g[0]; ++t)
      for (slice::const_iterator i = cb[k][t].begin();
           i != cb[k][t].end(); ++i)
        ASSERT_TRUE(i->empty());
  
}
TEST(Geometry, Periodicity){
  geometry::Geometry<DIM>::extents_t e;
  std::fill(e.begin(), e.end(), SIZE);
  geometry::Geometry<DIM> g(e);
  geometry::Geometry<DIM>::raw_pt_t n;
  std::fill(n.begin(), n.end(), 0);
  pt::Point<DIM> x_0 = g.mk_point(n);
  pt::Point<DIM> x = g.mk_point(n);
  for (; n[0] < SIZE; ++n[0], x += pt::Direction<DIM>(0))
    for (; n[1] < SIZE; ++n[1], x += pt::Direction<DIM>(1))
     for (; n[2] < SIZE; ++n[2], x += pt::Direction<DIM>(2))
       for (; n[3] < SIZE; ++n[3], x += pt::Direction<DIM>(3)) ;
  ASSERT_EQ(x, x_0);
}

/// Check if the reverse label lookup works
TEST(Geometry, ReverseLabelFn){
  geometry::Geometry<DIM>::extents_t e;
  std::fill(e.begin(), e.end(), SIZE);
  // make a funny looking lattice
  e[0] += 1;
  e[2] += 2;
  e[3] += 5;
  geometry::Geometry<DIM> g(e);
  geometry::Geometry<DIM>::raw_pt_t n, m;
  std::fill(n.begin(), n.end(), 0);
  for (; n[0] < SIZE; ++n[0])
    for (; n[1] < SIZE; ++n[1])
      for (; n[2] < SIZE; ++n[2])
        for (; n[3] < SIZE; ++n[3])
          ASSERT_EQ(n, g.coords(g.mk_point(n)));
}
