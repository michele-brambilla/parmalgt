#ifndef LOCAL_FIELD_HPP
#define LOCAL_FIELD_HPP

#include <vector>
#include <newQCDpt.h>
#include <Geometry.hpp>
#include <Point.hpp>
#include <Types.h>
#include <iostream>
#include <numeric>
#include <vector>
#ifdef MPI
#include <mpi.h>
#endif
#include <algorithm>
#ifdef _OPENMP
#include <omp.h>
#else
namespace fields {
  namespace detail {
    int omp_get_max_threads() { return 1; }
    int omp_get_thread_num() { return 0; }
  }
}
#endif


//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
///
///  Updated (MPI-enabled) versions of the perturbative fields.
///
///  The MPI enabled versions of the field variables have to take
///  care of communication.
///
///  \author Dirk Hesse <herr.dirk.hesse@gmail.com>
///  \date Mon Mar 26 14:49:53 2012
namespace fields {

  namespace detail {
    
    template <class Field_t> class inplace_add {
    public:
      static const int dim = Field_t::dim;
      inplace_add(const Field_t &F) : other(&F) { }
      void operator()(Field_t& F, const pt::Point<dim>& n) const {
        F[n] += (*other)[n];
      }
      static const int n_cb = 0;
    private:
      inplace_add(const inplace_add&) { } // no copy allowed
      Field_t const * other;
    };

   template <class Field_t> class inplace_sub {
    public:
      static const int dim = Field_t::dim;
      inplace_sub(const Field_t &F) : other(&F) { }
      void operator()(Field_t& F, const pt::Point<dim>& n) const {
        F[n] -= (*other)[n];
      }
     static const int n_cb = 0;
    private:
      inplace_sub(const inplace_sub&) { } // no copy allowed
      Field_t const * other;
    };

    template <class Field_t, class Scalar_t> 
    class inplace_smul {
    public:
      static const int dim = Field_t::dim;
      inplace_smul(const Scalar_t &s) : other(s) { }
      void operator()(Field_t& F, const pt::Point<dim>& n) const {
        F[n] *= other;
      }
      static const int n_cb = 0;
    private:
      inplace_smul(const inplace_smul&) { } // no copy allowed
      Scalar_t other;
    };


    //////////////////////////////////////////////////////
    //////////////////////////////////////////////////////
    ///  
    ///  The name makes sense for vectors: in this case is
    ///  "inner product" in strict sense. For matrices
    ///  it means whatever is the definition of operator*
    ///  
    ///  \author Michele Brambilla <mib.mic@gmail.com>
    ///  \date Thu Jan 17 10:40:05 2013
    template <class Field_t> class inner_prod {
    public:
      static const int dim = Field_t::dim;
      static const int n_cb = 0;
      inner_prod(const Field_t &F) : 
        result(omp_get_max_threads(), Cplx(0,0)), other(&F) { }

      void operator()(const Field_t& F, const pt::Point<dim>& n) {
        result[omp_get_thread_num()] += F[n] * (*other)[n];
      }

      Cplx reduce() const {
        return std::accumulate(result.begin(), result.end(), Cplx(0,0));
      }
    private:
      inner_prod(const inner_prod&) { } // no copy allowed
      std::vector<Cplx> result;
      Field_t const * other;
    };


    //////////////////////////////////////////////////////
    //////////////////////////////////////////////////////
    ///  
    ///  In case of vectors this the element-by-element
    ///  product (defined according to opertator^).
    ///  For matrices it means whatever is the definition 
    ///  of operator^
    ///  
    ///  \author Michele Brambilla <mib.mic@gmail.com>
    ///  \date Thu Jan 17 10:40:05 2013
    template <class Field_t> class prod {
    public:
      static const int dim = Field_t::dim;
      static const int n_cb = 0;
      prod(const Field_t &F) : 
        result(omp_get_max_threads(), Cplx(0,0)), other(&F) { }

      void operator()(const Field_t& F, const pt::Point<dim>& n) {
        result[omp_get_thread_num()] += F[n] ^ (*other)[n];
      }

      Cplx reduce() const {
        return std::accumulate(result.begin(), result.end(), Cplx(0,0));
      }
    private:
      prod(const prod&) { } // no copy allowed
      std::vector<Cplx> result;
      Field_t const * other;
    };


  }

  /// \defgroup MPI

  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  ///
  ///  Local gluon field.
  ///
  ///  By local we mean that it resides on one node. One
  ///  LocalGluonField instance is always associated with one MPI
  ///  process, but will generally spawn multiple threads such that
  ///  the advantages of a cluster using nodes with 
  ///  modern CPUs with several cores each may be fully exploited.
  ///
  ///  \tparam BGF The background field class to be used.
  ///  \tparam DIM The number of space-time dimensions
  ///  \tparam ORD The perturbative order.
  ///
  ///  \ingroup MPI
  ///
  ///  \author Dirk Hesse <herr.dirk.hesse@gmail.com>
  ///  \date Mon Mar 26 14:44:03 2012

  template <class F, int DIM>
  class LocalField {
  public:
    typedef F data_t;
    static const int dim = DIM;
    typedef std::vector< data_t > rep_t;
    typedef typename rep_t::iterator iterator;
    typedef typename rep_t::const_iterator const_iterator;
    typedef pt::Point<DIM> Point;
    typedef typename array_t<std::pair<int,int>, DIM>::Type neighbors_t;
    typedef typename geometry::Geometry<DIM>::extents_t extents_t;
    typedef typename geometry::Geometry<DIM>::raw_pt_t raw_pt;
    LocalField (const extents_t& e,
                const int& number_of_threads,
                const int& mpi_process_id,
                const neighbors_t& mpi_neighbors) : 
      g(e), rep(g.vol()), n_th(number_of_threads), pid(mpi_process_id),
      neighbors(mpi_neighbors){
#ifdef MPI
      /// constuct the buffers for communication
      for (int i = 0; i < DIM; ++i){
        send_buffer.push_back
          (std::make_pair ( std::vector<double> 
                            ( g.bnd_vol(i)*data_t::storage_size ),
                            std::vector<double> 
                            ( g.bnd_vol(i)*data_t::storage_size ) ) );
        rec_buffer.push_back
          (std::make_pair ( std::vector<double> 
                            ( g.bnd_vol(i)*data_t::storage_size ),
                            std::vector<double> 
                            ( g.bnd_vol(i)*data_t::storage_size ) ) );

      }
#endif
    }
    int vol() const { return g.vol(); }
    int extent(const int& i) const {
      return g[i];
    }
    
    data_t& operator[](const pt::Point<DIM> &n){
      return n.template deref<rep_t>(rep);
    }
    const data_t& operator[](const pt::Point<DIM> &n) const {
      return n.template deref<rep_t>(rep);
    }
    
    iterator begin() { return rep.begin(); };
    iterator end() { return rep.end(); };
    const_iterator begin() const { return rep.begin(); };
    const_iterator end() const { return rep.end(); };

    pt::Point<DIM> mk_point(const raw_pt& n){
      return g.mk_point(n);
    }
    geometry::SliceIterator<DIM, 0> mk_slice_iterator 
    (const pt::Direction<DIM> mu, const int& xi){
      return g.template mk_slice_iterator<0>(mu, xi);
    }
    geometry::SliceIterator<DIM, 2> mk_bulk_slice_iterator
    (const pt::Direction<DIM> mu, const int& xi){
      return g.mk_slice_iterator<2>(mu, xi, 1);
    }
    template <class M>
    void measure_on_slice(M& f, 
                          const pt::Direction<DIM>& d,
                          const int& xi){
      geometry::SliceIterator<DIM, 1> iter =
        g.template mk_slice_iterator<1>(d, xi, 0);
      while (iter.is_good()) f(*this, iter.yield());
    }
    template <class M>
    void measure_on_slice_with_bnd(M& f, 
                          const pt::Direction<DIM>& d,
                          const int& xi){
      geometry::SliceIterator<DIM, 0> iter =
        g.template mk_slice_iterator<0>(d, xi, 0);
      while (iter.is_good()) f(*this, iter.yield());
    }
    template <class M>
    void apply_on_slice_with_bnd(M& f, 
                          const pt::Direction<DIM>& d,
                          const int& xi){
      geometry::SliceIterator<DIM, 0> iter =
        g.template mk_slice_iterator<0>(d, xi, 0);
      while (iter.is_good()) f(*this, iter.yield());
    }

    template <class M>
    M& apply_on_timeslice(M& f, const int& t){
      // parallelize with a simple checker-board scheme ...
      typedef typename geometry::CheckerBoard<DIM, M::n_cb>::v_slice slice;
      typedef typename geometry::CheckerBoard<DIM, M::n_cb>::v_bin bin;
      geometry::CheckerBoard<DIM, M::n_cb> cb(g);
      for (typename slice::const_iterator s = cb[t].begin();
           s != cb[t].end(); ++s){
        // here, we have to do a nasty workaround
        // parallel for does not like lists, which is what I used in
        // the CheckerBoard class. This is e.g. great for the test
        // case of the CB class and foremost for the case that we have
        // some overlap at the edges (i.e. if CHECKER_BOARD_SIZE % L
        // != 0). Here, however a vector would be nice. Let's see...
        // FIXME: This is nasty!!
        //std::vector<pt::Point<DIM> > v(s->begin(), s->end());
        int N = s->size();
#pragma omp parallel for
        for (int i = 0; i < N; ++i)
          f(*this, (*s)[i]);
      }
      return f;
    }
    template <class M>
    M& apply_on_timeslice(M& f, const int& t) const {
      // parallelize with a simple checker-board scheme ...
      typedef typename geometry::CheckerBoard<DIM, M::n_cb>::v_slice slice;
      typedef typename geometry::CheckerBoard<DIM, M::n_cb>::v_bin bin;
      static geometry::CheckerBoard<DIM, M::n_cb> cb(g);
      for (typename slice::const_iterator s = cb[t].begin();
           s != cb[t].end(); ++s){
        int N = s->size();
#pragma omp parallel for
        for (int i = 0; i < N; ++i)
          f(*this, (*s)[i]);
      }
      return f;
    }
    template <class M>
    M& apply_everywhere(M& f){
      for (int t = 0; t < g[0]; ++t)
        apply_on_timeslice(f, t);
      return f;
    }
    template <class M>
    M& apply_everywhere(M& f) const {
      for (int t = 0; t < g[0]; ++t)
        apply_on_timeslice(f, t);
      return f;
    }
    template <class M>
    const M& apply_everywhere(const M& f){
      for (int t = 0; t < g[0]; ++t)
        apply_on_timeslice(f, t);
      return f;
    }
    template <class M>
    void apply_everywhere_serial(M& f){
      // doesn't work because bind1st makes a copy
      // we don't want that because of local data members
      // such as the fstreams in the disk reader/writer
      //std::for_each(g.begin(), g.end(), std::bind1st(f, *this));
      for(pt::Point<DIM> n = g.begin(), e = g.end(); n != e; ++n)
        f(*this, n);
    }
#ifdef MPI
    MPI_Request test_send_fwd_z(){
      write_slice_to_buffer(pt::Direction<DIM>(3), 4,
                            send_buffer[3].second);
      MPI_Request r;
      MPI_Isend(&send_buffer[3].second[0],
                send_buffer[3].second.size(),
                MPI_DOUBLE, neighbors[3].second, 0, MPI_COMM_WORLD, &r);
      return r;
    }
    void test_rec_bkw_z(){
      MPI_Status status;
      MPI_Recv(&rec_buffer[3].first[0],
      rec_buffer[3].first.size(),
               MPI_DOUBLE, neighbors[3].first, 0, MPI_COMM_WORLD,
               &status);
      read_slice_from_buffer(pt::Direction<DIM>(3), 0,
                             rec_buffer[3].first);
    }
#endif

    // Added: MB on Fri 2, 2012

    LocalField& operator+=(LocalField& other) {
      apply_everywhere(detail::inplace_add<LocalField>(other));
      return *this;
    }
    LocalField& operator-=(LocalField& other) {
      apply_everywhere(detail::inplace_sub<LocalField>(other));
      return *this;
    }
    LocalField& operator+=(const LocalField& other) {
      apply_everywhere(detail::inplace_add<LocalField>(other));
      return *this;
    }
    LocalField& operator-=(const LocalField& other) {
      apply_everywhere(detail::inplace_sub<LocalField>(other));
      return *this;
    }

    LocalField operator+(const LocalField& other) {
      LocalField result(*this);
      return result+=other;
    }
    LocalField operator-(const LocalField& other) {
      LocalField result(*this);
      return result-=other;
    }

    template<class C>
    LocalField& operator*=(C& other) {
      apply_everywhere(detail::inplace_smul<LocalField, C>(other));
      return *this;
    }
    template<class C>
    LocalField operator*(const C& other) const {
      LocalField result(*this);
      return (result *= other);
    }
    // DH: maybe call this prod instead?
    Cplx operator*(const LocalField& other) const { 
      detail::inner_prod<LocalField> p(other);
      return apply_everywhere(p).reduce();
    }

    ///  \author Michele Brambilla <mib.mic@gmail.com>
    ///  \date Thu Jan 10 20:47:14 2013
    Cplx operator^(const LocalField& other) const { 
      detail::prod<LocalField> p(other);
      return apply_everywhere(p).reduce();
    }
    
  private:
    geometry::Geometry<DIM> g;
    rep_t rep;
    int n_th; // number of threads
    int pid; // MPI process id
    neighbors_t neighbors;
    /// Buffer for communication
    std::vector< std::pair<std::vector<double>, 
                           std::vector<double> > > send_buffer;
    std::vector< std::pair<std::vector<double>, 
                           std::vector<double> > > rec_buffer;
    void write_slice_to_buffer(const pt::Direction<DIM>& mu, const int &xi,
                               std::vector<double>& buff){
      std::vector<double>::iterator i = buff.begin();
      geometry::SliceIterator<DIM, 0> iter = 
        g.template mk_slice_iterator<0>(mu, xi);
      while (iter.is_good()) rep[iter.yield()].buffer(i);
    }
    void read_slice_from_buffer(const pt::Direction<DIM>& mu, const int &xi,
                                std::vector<double>& buff){
      std::vector<double>::const_iterator i = buff.begin();
      geometry::SliceIterator<DIM, 0> iter = 
        g.template mk_slice_iterator<0>(mu, xi);
      while (iter.is_good()) rep[iter.yield()].unbuffer(i);
    }
  };
}

#endif
