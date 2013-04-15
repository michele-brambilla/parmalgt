#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <Types.h>
#include <Point.hpp>
#include <math.h>
#include <list>

namespace geometry {

  namespace detail {
    template <int N>
    inline int product (const typename array_t<int, N>::Type &v){
      int p = 1;
      for (typename array_t<int, N>::Type::const_iterator i = v.begin(); i !=
             v.end(); ++i) p *= *i;
      return p;
    }    
  }

  template<int DIM, int SKIP> class SliceIterator;

  
  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  ///
  ///  Lattice geometry helper class.
  ///
  ///  This class stores
  ///     - A list of neighbors (Geometry::neighbors).
  ///     - The extents of the lattice (Geometry::extents).
  ///     - The DIM-volume Geometry::V.
  ///     - The volumes of the boundaries (Geometry::bnd_vols).
  ///
  ///  \tparam DIM Space-time dimensions.
  ///
  ///  \ingroup MPI
  ///  \author Dirk Hesse <herr.dirk.hesse@gmail.com>
  ///  \date Mon Mar 26 12:12:02 2012
  template <int DIM>
  class Geometry {
  public:
    /// Type for integer DIM-Vectors.
    /// Used for the construction of a point n = (n1, n2, ...).
    /// For the day-to day use, you could write a wrapper function
    /// or consider the begin() and end() methods to swipe the whole
    /// lattice .
    typedef typename array_t<int, DIM>::Type raw_pt_t;
    /// Type to pass the lattice extents/
    /// The same as raw_pt_t, an integer DIM-vector.
    typedef typename array_t<int, DIM>::Type extents_t;
    /// Type to store the volumes of the boundaries.
    typedef extents_t vols_t;
    /// Constructor.
    /// Only needs lattice extents as input. We recommend using a
    /// helper function to create the parameter.
    Geometry (const extents_t &lattice_extents)
    : extents(lattice_extents),
      V(detail::product<DIM>(lattice_extents)),
      neighbors (detail::product<DIM>(lattice_extents), 
                  typename array_t<int, 2*DIM>::Type()) { 
      // initialize neighbors
      fill_neighbors();
      // boundary volumes
      for (int i = 0; i < DIM; ++i)
        bnd_vols[i] = V / extents[i];
    }
    extents_t get_extents() const { return extents; }
    /// Generate an iterator over a given slice.
    /// The SliceIterator returned will ierate all lattice points
    /// \f$n\f$, keeping \f$n_\mu = x_i\f$ fixed.
    /// For all other directions, the iteration range is
    /// \f[ 0 \leq x_\nu < L_\nu, \quad \nu = 0,\ldots, DIM - 1, \quad  \nu \neq \mu\f]
    ///
    ///
    /// \param mu Called \f$\mu\f$ above.
    /// \param xi Called \f$n_i\f$ above.
    template <int N>
    SliceIterator<DIM, N> mk_slice_iterator (const pt::Direction<DIM> mu,
                                             const int& xi,
                                             const int& x_start = 0) const {
      raw_pt_t n;
      std::fill(n.begin(), n.end(), x_start);
      n[mu] = xi;
      return SliceIterator<DIM, N>(mk_point(n), mu, extents);
    }
    
    /// Generate an iterator over a given slice, omitting the boundaries.
    /// The SliceIterator returned will ierate all lattice points
    /// \f$n\f$, keeping \f$n_\mu = x_i\f$ fixed.
    /// For all other directions, the iteration range is
    /// \f[ 1 \leq x_\nu < L_\nu - 1, \quad \nu = 0,\ldots, DIM - 1, \quad  \nu \neq \mu\f]
    ///
    ///
    /// \param mu Called \f$\mu\f$ above.
    /// \param xi Called \f$n_i\f$ above.
    template <int N>
    SliceIterator<DIM, 2> mk_vol_iterator (const pt::Direction<DIM> mu,
                                          const int& xi){
      raw_pt_t n;
      std::fill(n.begin(), n.end(), 1);
      n[mu] = xi;
      
    }
    pt::Point<DIM> mk_point(const raw_pt_t &n) const {
      return pt::Point<DIM>(mk_label(n), neighbors.begin());
    }
    raw_pt_t coords(const pt::Point<DIM>& p) const {
      return reverse_label(p);
    }
    pt::Point<DIM> begin() const {
      return pt::Point<DIM>(0, neighbors.begin());
    }
    pt::Point<DIM> end() const {
      return pt::Point<DIM>(V, neighbors.begin());
    }
    int bnd_vol(const int& dim) const {
      return bnd_vols[dim];
    }
    int vol() const { return V; }

    //////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////
    ///
    ///  Checker board scheme bin numbers.
    ///
    ///  We implement a checker board scheme with cube size N. It
    ///  consists of two steps.
    ///  1) Divde the lattice in hypercubes of size N^D. Label these
    ///     alternating (0) and (1).
    ///  2) Within the cubes, number your elements consistently, and
    ///     combine this witht he result from 1) into the final bin
    ///     number. Thus, a point in a cube labeled with (0) and with
    ///     number 8 within this hypercube will end up in bin (0)->8.
    ///  Thus, in total we have 2*N^D bins. In the case N=1 this
    ///  reduces to the usual checkerboard scheme.
    ///
    ///  added by DH, 31st Oct. 2012
    template <int N>
    int bin(const pt::Point<DIM>& p) const {
      raw_pt_t n = coords(p);
      // A) is p in a black (0) or a white (1) block?
      int A = 0;
      for (int i = 0; i < DIM; ++i)
        A = (A + (int(floor(n[i]/N)) % 2)) % 2;
      // B) which position within the block does p have?
      int vol = 1, B = 0;
      for (int i = 0; i < DIM; ++i, vol *= N)
        B += vol*(n[i] % N);
      return A*pow(N, DIM) + B;
    }


    /// get the total number of bins with a checkerboard scheme of
    /// cube size N
    template <int N>
    int n_bins() const{
      return 2*pow(N, DIM);
    }

    //////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////
    ///
    ///  Get the extend in a certain direction.
    ///
    ///  \author Dirk Hesse <herr.dirk.hesse@gmail.com>
    ///  \date Mon Jun 11 18:33:23 2012
    int operator[](const int& i) const { return extents[i]; }
  private:
    /// Local lattice extends.
    extents_t extents;
    /// Stores the nearest neighbors.  Given a point \f$n\f$,
    /// neighbors[n][DIM + mu] is the coordinate of 
    /// \f$ n + \hat \mu \f$, neighbors[n][mu] the one of 
    /// \f$ n - \hat \mu \f$, with \f$ \mu = 0,\ldots, DIM\f$.

    std::vector<typename array_t<int, 2*DIM>::Type> neighbors;
    /// Voulumes of the boundaries
    vols_t bnd_vols;
    /// DIM-voulume of the local lattice
    int V;
    /// Label (integer value) belonging to a give point
    int mk_label(const raw_pt_t &n) const {
      int vol = 1, label = 0;
      for (int i = DIM-1; i >= 0; --i){
        label += ((n[i] + extents[i]) % extents[i]) * vol;
        vol *= extents[i];
      }
      return label;
    }
    /// Get the coordinates from a given label
    raw_pt_t reverse_label(const int& n) const {
      int vol = 1;
      raw_pt_t x;
      for (int i = DIM - 1; i >= 0; --i){
        x[i] = (n / vol) % extents[i];
        vol *= extents[i];
      }
      return x;
    }
    /// helper function to recursively fill neighbors
    void fill_neighbors (const int& n = DIM - 1, raw_pt_t x = raw_pt_t()){
      if (!n)
        for (int xi = 0; xi < extents[0]; ++xi){
          x[0] = xi;
          int l = mk_label(x);
          for(pt::Direction<DIM> mu; mu.is_good(); ++mu){
            x[mu]++;
            neighbors[l][DIM + (int)mu] = mk_label(x);
            x[mu] -= 2;
            neighbors[l][mu] = mk_label(x);
            x[mu]++;
          }
        }
      else{
        for (int xi = 0; xi < extents[n]; ++xi){
          x[n] = xi;
          fill_neighbors (n-1, x);
        }
      }
    }
  };
  

  namespace new_iter {
    ////////////////////////////////////////////////////////////
    //
    //  Increment fuctions.
    //
    //  These are policies how to increment the individual
    /// dimensions when iterating
    //
    //  \date      Tue Apr  9 16:27:53 2013
    //  \author    Dirk Hesse <dirk.hesse@fis.unipr.it>

    template <int DIM>
    struct ConstPolicy {
      ConstPolicy (int) {}
      void next(const pt::Point<DIM>&, const pt::Direction<DIM>&) { }
      void reset(void) { }
      bool overflow(void) { return true; }
    };

    template <int DIM>
    struct PeriodicPolicy {
      int curr, max;
      PeriodicPolicy (int m) : curr(0), max(m) { }
      void next(pt::Point<DIM> &p , const pt::Direction<DIM> &mu) {
        curr++;
        p += mu;
      }
      void reset(void) { curr = 0; }
      bool overflow(void) { return curr >= max; }
    };

    template <int DIM, class P, class N>
    struct Pair {
      P policy;
      N next;
      bool good;
      template <class InputIterator>
      Pair (InputIterator L) :
      policy(*(L++)), next(L), good(true) { }
      bool advance(pt::Point<DIM>& p,
                   pt::Direction<DIM> mu =
                   pt::Direction<DIM>(0)){
        policy.next(p, mu);
        if (policy.overflow()){
          policy.reset();
          good = next.advance(p, ++mu);
        }
        return good;
      }
      bool is_good() { return good; }
    };

    template <int DIM>
    struct End {
       template <class InputIterator>
      End(InputIterator) { }
      bool advance(pt::Point<DIM>&,
                   pt::Direction<DIM>&) { return false; }
      bool is_good() { return false; }
    };

#define RAW_TWO_POLICY_ITER(I, A, B)             \
    Pair<I, A<I>, Pair<I, B<I>, End<I> > >
#define TWO_POLICY_ITER(A, B)                   \
    RAW_TWO_POLICY_ITER(2, A, B)
#define THREE_POLICY_ITER(A, B, C)                \
    Pair< 3, A<3>, RAW_TWO_POLICY_ITER(3, B, C) >

    typedef TWO_POLICY_ITER(PeriodicPolicy,
                            PeriodicPolicy) PeriodicTwoDimIter;
    typedef THREE_POLICY_ITER(PeriodicPolicy,
                              PeriodicPolicy,
                              ConstPolicy) PeriodicConstThreeDimIter;
  }

  // specializations for four dimensions ...

  template <> template <>
  int Geometry<4>::bin<0>(const pt::Point<4>&) const { return 0; }
  template <> template <>
  int Geometry<4>::n_bins<0>() const { return 1; }

  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  ///
  ///  Class to facilitate the iteration over a slice of a Geometry.
  ///
  ///  Given a point \f$x\f$, and a direction \f$\mu\f$, yield will
  ///  successively return a sequence of points \f$x^{(i)}\f$, with
  ///  the \f$\mu\f$ component fixed \f$x^{(i)}_\mu = x_\mu\f$, but
  ///  all the other components sweeping the whole lattice extent.
  ///
  ///  \tparam DIM Number of space-time dimensions.
  ///
  ///  \ingroup MPI
  ///
  ///  \author Dirk Hesse <herr.dirk.hesse@gmail.com>
  ///  \date Mon Mar 26 12:58:50 2012

  template<int DIM, int SKIP>
  class SliceIterator {
  public:
    /// Constructor.
    /// \param x The starting point, called \f$x\f$ in the class
    ///          description.
    /// \param mu The direction to keep fixed, called \f$\mu\f$ int
    ///          the class escription.
    /// \param e Extents of the underlying Geometry.
    SliceIterator (const pt::Point<DIM> &x,
                   const pt::Direction<DIM> mu,
                   const typename Geometry<DIM>::extents_t& e,
                   const int& stepsize = 1) :
      x_current(x), mu_exclude(mu), extents(e), counters(),
      good_flag(true) { }

    /// Return the next point.
    /// Note you should not call this if is_good() returns false.
    ///
    /// Typical use case:
    /// 
    pt::Point<DIM> yield() {
      pt::Point<DIM> result = x_current;
      int mu = (0 == mu_exclude) ? 1 : 0;
      do {
          x_current += pt::Direction<DIM>(mu);
          counters[mu]++;
        if (counters[mu] == extents[mu] - SKIP){
          counters[mu] = 0;
          // this uses periodicity !!
          for (int n = 0; n < SKIP; ++n)
            x_current += pt::Direction<DIM>(mu);
          if (++mu == mu_exclude) ++mu;
        }
        else break;
      } while (mu <= DIM);
      if (mu == DIM) good_flag = false;
      return result;
    }

    bool is_good() const {
      return good_flag;
    }
  private:
    pt::Point<DIM> x_current;
    pt::Direction<DIM> mu_exclude;
    typename Geometry<DIM>::extents_t extents;
    typename array_t<int, DIM>::Type counters;
    bool good_flag;
    int step;
  };

  /////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////
  ///
  ///  Checkerboard scheme for safe parallel application of operators
  ///  etc. For the moment, I gave the class two copies of the lists
  ///  of bins, one using c++ list, one using vectors. The list one is
  ///  mainly for ease of construction and for use in the unit
  ///  tests. The other one is to speed up the iteration at run-time
  ///  in LocalField.
  ///
  ///  DH, 7th Nov. 2012
  template <int D, int N>
  class CheckerBoard {
  public:
    typedef pt::Point<D> Point; // single point
    // first, we use a list for fast generation of the bins
    typedef std::list<Point> l_bin; // one bin
    typedef std::vector<l_bin> l_slice; // bins in a time slice
    typedef std::vector<l_slice> l_lattice; // all time slices
    // in a second step, we replace that by vectors for later
    // convenience
    typedef std::vector<Point> v_bin; // one bin
    typedef std::vector<v_bin> v_slice; // bins in a time slice
    typedef std::vector<v_slice> v_lattice; // all time slices
    explicit CheckerBoard (const Geometry<D>& g) :
      lat(g[0] + 1, l_slice(g.template n_bins<N>())), v_lat(g[0] + 1) {
      for (int t = 0; t <= g[0]; ++t){
        geometry::SliceIterator<D, 0> iter =
          g.template mk_slice_iterator<0>(pt::Direction<D>(0), t, 0);
        while (iter.is_good()){
          Point p = iter.yield();
          lat[t].at(g.template bin<N>(p)).push_back(p);
        }
        for (typename l_slice::const_iterator i = lat[t].begin();
             i != lat[t].end(); ++i)
          v_lat[t].push_back(v_bin(i->begin(), i->end()));
      }
    }
    const v_slice& operator[](const int& i) const {
      return v_lat[i];
    }
    // as operator[] just for non-const access. this is for testing
    // purposes only and should not be used 'in the wild' use the
    // constant operator above instead.
    l_slice& nca(const int& i){
      return lat[i];
    }
  private:
    l_lattice lat;
    v_lattice v_lat;
  };
}

#endif
