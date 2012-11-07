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
            ASSERT_NE(g.bin<1>(x), g.bin<1>(x + mu));
            ASSERT_EQ(g.bin<1>(x), g.bin<1>(x + mu + mu));
            for (pt::Direction<DIM> nu; nu.is_good(); ++nu){
              ASSERT_NE(g.bin<2>(x), g.bin<2>(x + mu + nu));
              for (pt::Direction<DIM> rho; rho.is_good(); ++rho)
                ASSERT_NE(g.bin<2>(x), g.bin<2>(x + mu + nu + rho));
            }
            ASSERT_EQ(g.bin<2>(x), g.bin<2>(x + mu + mu + mu + mu));
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
  geometry::CheckerBoard<DIM,0> cb0(g);
  geometry::CheckerBoard<DIM,1> cb1(g);
  geometry::CheckerBoard<DIM,2> cb2(g);
  geometry::CheckerBoard<DIM,3> cb3(g);
  typedef geometry::CheckerBoard<DIM,0>::slice slice0;
  typedef geometry::CheckerBoard<DIM,1>::slice slice1;
  typedef geometry::CheckerBoard<DIM,2>::slice slice2;
  typedef geometry::CheckerBoard<DIM,3>::slice slice3;
  typedef geometry::CheckerBoard<DIM,0>::bin bin0;
  typedef geometry::CheckerBoard<DIM,1>::bin bin1;
  typedef geometry::CheckerBoard<DIM,2>::bin bin2;
  typedef geometry::CheckerBoard<DIM,3>::bin bin3;
  geometry::Geometry<DIM>::raw_pt_t n;
  std::fill(n.begin(), n.end(), 0);
  for (n[0] = 0; n[0] < e[0]; ++n[0])
    for (n[1] = 0; n[1] < e[1]; ++n[1])
      for (n[2] = 0; n[2] < e[2]; ++n[2])
        for (n[3] = 0; n[3] < e[3]; ++n[3]){
          pt::Point<DIM> x = g.mk_point(n);
          ///////////////////////////////
          /// N = 0
          bool found = false;
          // note: use nca instead of operator[] below for
          // non-constant access which we will use to erase the found
          // elements to assure we don't have duplicates or somehting
          // the like (which would be very strange!)
          for (slice0::iterator i = cb0.nca(n[0]).begin();
                 i != cb0.nca(n[0]).end() && !found; ++i){
            bin0::iterator f = std::find(i->begin(), i->end(), x);
            if (f != i->end()){
              found = true;
              i->erase(f);
            }
          }
          ASSERT_EQ(found, true);
          ///////////////////////////////
          /// N = 1
          found = false;
          // note: use nca instead of operator[] below for
          // non-constant access which we will use to erase the found
          // elements to assure we don't have duplicates or somehting
          // the like (which would be very strange!)
          for (slice1::iterator i = cb1.nca(n[0]).begin();
                 i != cb1.nca(n[0]).end() && !found; ++i){
            bin1::iterator f = std::find(i->begin(), i->end(), x);
            if (f != i->end()){
              found = true;
              i->erase(f);
            }
          }
          ASSERT_EQ(found, true);
          ///////////////////////////////
          /// N = 2
          found = false;
          // note: use nca instead of operator[] below for
          // non-constant access which we will use to erase the found
          // elements to assure we don't have duplicates or somehting
          // the like (which would be very strange!)
          for (slice2::iterator i = cb2.nca(n[0]).begin();
                 i != cb2.nca(n[0]).end() && !found; ++i){
            bin2::iterator f = std::find(i->begin(), i->end(), x);
            if (f != i->end()){
              found = true;
              i->erase(f);
            }
          }
          ASSERT_EQ(found, true);
          ///////////////////////////////
          /// N = 3
          found = false;
          // note: use nca instead of operator[] below for
          // non-constant access which we will use to erase the found
          // elements to assure we don't have duplicates or somehting
          // the like (which would be very strange!)
          for (slice3::iterator i = cb3.nca(n[0]).begin();
                 i != cb3.nca(n[0]).end() && !found; ++i){
            bin3::iterator f = std::find(i->begin(), i->end(), x);
            if (f != i->end()){
              found = true;
              i->erase(f);
            }
          }
          ASSERT_EQ(found, true);
        }
  // check if we hit all points ..
  ///////////////////////////////
  /// N = 0
  for (int t = 0; t < g[0]; ++t)
    for (slice0::const_iterator i = cb0[t].begin();
         i != cb0[t].end(); ++i)
      ASSERT_TRUE(i->empty());
  ///////////////////////////////
  /// N = 1
  for (int t = 0; t < g[0]; ++t)
    for (slice1::const_iterator i = cb1[t].begin();
         i != cb1[t].end(); ++i)
      ASSERT_TRUE(i->empty());
  ///////////////////////////////
  /// N = 2
  for (int t = 0; t < g[0]; ++t)
    for (slice2::const_iterator i = cb2[t].begin();
         i != cb2[t].end(); ++i)
      ASSERT_TRUE(i->empty());
  ///////////////////////////////
  /// N = 3
  for (int t = 0; t < g[0]; ++t)
    for (slice3::const_iterator i = cb3[t].begin();
         i != cb3[t].end(); ++i)
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
