#include "Background.h"
#include "gtest/gtest.h"
#include "Point.hpp"
#include "Geometry.hpp"

const int DIM = 4;
const int SIZE = 4;


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
	    ASSERT_EQ(g.mk_point(x), g.mk_point(n) - (-mu));
            ASSERT_NE(g.mk_point(x), g.mk_point(n));
            x[mu] -= 2;
            ASSERT_EQ(g.mk_point(x), g.mk_point(n) - mu);
            ASSERT_EQ(g.mk_point(x), g.mk_point(n)  + (-mu));
            ASSERT_NE(g.mk_point(x), g.mk_point(n));
        }
};

TEST(SliceIteratorTest, CrossCheckKnownValues){
  geometry::Geometry<DIM>::extents_t e;
  std::fill(e.begin(), e.end(), SIZE);
  geometry::Geometry<DIM> g(e);
  geometry::Geometry<DIM>::raw_pt_t n;
  // create a list with all points at time SIZE / 2
  typedef std::list<pt::Point<DIM> > list_t;
  list_t known;
  n[0] = SIZE / 2;
  for (n[1] = 0; n[1] < SIZE; ++n[1])
    for (n[2] = 0; n[2] < SIZE; ++n[2])
      for (n[3] = 0; n[3] < SIZE; ++n[3])
        known.push_back(g.mk_point(n));
  // try to re-create that list using a SliceIterator...
  geometry::new_iter::TimeSliceIter s_iter(e.begin());
  pt::Point<DIM> x = g.mk_point(n);
  do {
    // is iter.yield() in the list?
    //pt::Point<DIM> tmp = s_iter.yield();
    list_t::iterator findr = std::find(known.begin(),
                                       known.end(),
                                       x);
    geometry::Geometry<DIM>::raw_pt_t m(g.coords(x));
    ASSERT_TRUE(findr != known.end()) << m[0] << ","
                                      << m[1] << ","
                                      << m[2] << ","
                                      << m[3];
    known.erase(findr);
    s_iter.advance(x);
  } while (s_iter.is_good());
  ASSERT_EQ(0, known.size());
}

TEST(SliceIteratorTest, CrossCheckKnownValuesBulkOnly){
  geometry::Geometry<DIM>::extents_t e;
  std::fill(e.begin(), e.end(), SIZE);
  geometry::Geometry<DIM> g(e);
  geometry::Geometry<DIM>::raw_pt_t n;
  // create a list with all points at time SIZE / 2
  typedef std::list<pt::Point<DIM> > list_t;
  list_t known;
  n[0] = SIZE / 2;
  for (n[1] = 0; n[1] < SIZE; ++n[1])
    for (n[2] = 0; n[2] < SIZE; ++n[2])
      for (n[3] = 1; n[3] < SIZE-1; ++n[3])
        known.push_back(g.mk_point(n));
  n[1] = 1; n[2] = 1; n[3] = 1;
  // try to re-create that list using a SliceIterator...
  //geometry::SliceIterator<DIM, 2>
  //  s_iter(g.mk_point(n), pt::Direction<DIM>(0), e);
  geometry::new_iter::ZBulkTimeSliceIter s_iter(e.begin());
  pt::Point<DIM> tmp = g.mk_point(n);
  do {
    // is iter.yield() in the list?
    list_t::iterator findr = std::find(known.begin(),
                                       known.end(),
                                       tmp);
    geometry::Geometry<DIM>::raw_pt_t m(g.coords(tmp));
    ASSERT_TRUE(findr != known.end()) << m[0] << ","
                                      << m[1] << ","
                                      << m[2] << ","
                                      << m[3];

    known.erase(findr);
    s_iter.advance(tmp);
  } while (s_iter.is_good());
  ASSERT_EQ(0, known.size());
}

TEST(NewIterator, TwoPeriodicPol){
  geometry::Geometry<2>::extents_t e, x;
  e[1] = 3; e[0] = 3;
  x[0] = 0; x[1] = 0;
  geometry::new_iter::PeriodicTwoDimIter i(e.begin());
  geometry::Geometry<2> g(e);
  pt::Point<2> p = g.mk_point(x);
  int count = 0;
  do {
    ASSERT_EQ(g.coords(p)[0], count / 3);
    ASSERT_EQ(g.coords(p)[1], count++ % 3);
    i.advance(p);
  } while (i.is_good());
}

TEST(NewIterator, ThreePeriodicConstPol){
  geometry::Geometry<3>::extents_t e, x;
  e[1] = 3; e[0] = 3; e[2] = 3;
  x[0] = 2; x[1] = 0; x[2] = 0;
  geometry::new_iter::PeriodicConstThreeDimIter i(e.begin());
  geometry::Geometry<3> g(e);
  pt::Point<3> p = g.mk_point(x);
  int count = 0;
  do {
    ASSERT_EQ(g.coords(p)[0], 2);
    ASSERT_EQ(g.coords(p)[1], count / 3);
    ASSERT_EQ(g.coords(p)[2], count++ % 3);
    i.advance(p);
  } while (i.is_good());
}

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
  typedef geometry::CheckerBoard<DIM,0>::l_slice slice0;
  typedef geometry::CheckerBoard<DIM,1>::l_slice slice1;
  typedef geometry::CheckerBoard<DIM,2>::l_slice slice2;
  typedef geometry::CheckerBoard<DIM,3>::l_slice slice3;
  typedef geometry::CheckerBoard<DIM,0>::l_bin bin0;
  typedef geometry::CheckerBoard<DIM,1>::l_bin bin1;
  typedef geometry::CheckerBoard<DIM,2>::l_bin bin2;
  typedef geometry::CheckerBoard<DIM,3>::l_bin bin3;
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
    for (slice0::iterator i = cb0.nca(t).begin();
         i != cb0.nca(t).end(); ++i)
      ASSERT_TRUE(i->empty());
  ///////////////////////////////
  /// N = 1
  for (int t = 0; t < g[0]; ++t)
    for (slice1::iterator i = cb1.nca(t).begin();
         i != cb1.nca(t).end(); ++i)
      ASSERT_TRUE(i->empty());
  ///////////////////////////////
  /// N = 2
  for (int t = 0; t < g[0]; ++t)
    for (slice2::iterator i = cb2.nca(t).begin();
         i != cb2.nca(t).end(); ++i)
      ASSERT_TRUE(i->empty());
  ///////////////////////////////
  /// N = 3
  for (int t = 0; t < g[0]; ++t)
    for (slice3::iterator i = cb3.nca(t).begin();
         i != cb3.nca(t).end(); ++i)
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
