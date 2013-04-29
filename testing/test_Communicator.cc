#include <gtest/gtest.h>
#include <LocalField.hpp>
#include <cstdlib>
#include <Helper.h>
#include <Kernels.hpp>

template <class Field_t>
struct randomize {
  randomize(const int& seed) { std::srand(seed); }
  static const int dim = Field_t::dim;
  void operator()(Field_t& F, const pt::Point<dim>& n) const {
    F[n] = rand();
  }
  static const int n_cb = 0;
};

template <class Field_t>
struct cfld {
  cfld(const int& value) : val(value) { }
  static const int dim = Field_t::dim;
  void operator()(Field_t& F, const pt::Point<dim>& n) const {
    F[n] = val;
  }
  static const int n_cb = 0;
  int val;
};

// TEST(Comunicator, ApplyOnTimesliceIteratorTest){
//   typedef fields::LocalField<int, 4> intField;
//   typedef typename intField::raw_pt raw_pt_t;
//   typedef typename std::vector<int> vec_t;
//   intField::extents_t e;
//   intField::neighbors_t nn;
//   std::fill(e.begin(), e.end(), 4);
//   intField A(e,1,0,nn);
//   cfld<intField> c(A.get_communicator().rank()+1);
//   A.apply_on_timeslice(c,1);
//   vec_t v, vtest;
//   for(int z = 0; z < 4; ++z) {
//     v.push_back(A[A.mk_point(raw_pt_t{1,0,0,z})]);
//     vtest.push_back(0);
//   }

//   for( int rank = 0; rank < A.get_communicator().numprocs(); ++rank) {
//     for( vec_t::iterator it = vtest.begin()+1; it != vtest.end()-1; ++it) *it = rank+1;

//     if( A.get_communicator().rank() == rank)
//       for( vec_t::iterator it = vtest.begin(), iv = v.begin(); 
// 	   it != vtest.end(); ++it, ++iv)
// 	ASSERT_EQ(*it, *iv);
    
//     MPI_Barrier(MPI_COMM_WORLD);
//     std::cout << std::endl;
//   }
// }

// TEST(Comunicator, EverywhereIteratorTest){
//   typedef fields::LocalField<int, 4> intField;
//   typedef typename intField::raw_pt raw_pt_t;
//   typedef typename std::vector<int> vec_t;
//   intField::extents_t e;
//   intField::neighbors_t nn;
//   std::fill(e.begin(), e.end(), 4);
//   intField A(e,1,0,nn);
//   cfld<intField> c(A.get_communicator().rank()+1);
//   A.apply_everywhere(c);
//   vec_t v, vtest;
//   for(int z = 0; z < 4; ++z) {
//     v.push_back(A[A.mk_point(raw_pt_t{0,0,0,z})]);
//     vtest.push_back(0);
//   }

//   for( int rank = 0; rank < A.get_communicator().numprocs(); ++rank) {
//     for( vec_t::iterator it = vtest.begin()+1; it != vtest.end()-1; ++it) *it = rank+1;

//     if( A.get_communicator().rank() == rank)
//       for( vec_t::iterator it = vtest.begin(), iv = v.begin(); 
// 	   it != vtest.end(); ++it, ++iv)
// 	ASSERT_EQ(*it, *iv);
    
//     MPI_Barrier(MPI_COMM_WORLD);
//     std::cout << std::endl;
//   }
// }


TEST(Comunicator, Buffer){
  typedef fields::LocalField<int, 4> intField;
  typedef typename intField::raw_pt raw_pt_t;
  typedef typename std::vector<int> vec_t;
  intField::extents_t e;
  intField::neighbors_t nn;
  std::fill(e.begin(), e.end(), 4);
  intField A(e,1,0,nn);
  cfld<intField> c(A.get_communicator().rank()+1);
  A.apply_everywhere(c);
  A.apply_everywhere(c);

  ASSERT_EQ((A.get_communicator().send_buff()[0].first).size(), e[1]*e[2]*e[3] );
  ASSERT_EQ((A.get_communicator().send_buff()[3].first).size(), e[1]*e[2] );
  
  for( int rank = 0; rank < A.get_communicator().numprocs(); ++rank) {
    // test send buffer    
    if( A.get_communicator().rank() == rank) {
      for( std::vector<int>::const_iterator it = A.get_communicator().send_buff()[3].first.begin();
	   it != A.get_communicator().send_buff()[3].first.end(); ++it)
	ASSERT_EQ(*it, (A.get_communicator().rank()+1));
      for( std::vector<int>::const_iterator it = A.get_communicator().send_buff()[3].second.begin();
	   it != A.get_communicator().send_buff()[3].second.end(); ++it)
	ASSERT_EQ(*it, (A.get_communicator().rank()+1));
    }

    MPI_Barrier(MPI_COMM_WORLD);
  }
}



TEST(Comunicator, SendRecv){
  typedef fields::LocalField<int, 4> intField;
  typedef typename intField::raw_pt raw_pt_t;
  typedef typename std::vector<int> vec_t;
  intField::extents_t e;
  intField::neighbors_t nn;
  std::fill(e.begin(), e.end(), 4);
  intField A(e,1,0,nn);
  cfld<intField> c(A.get_communicator().rank()+1);

  const int numprocs = A.get_communicator().numprocs();
  typedef typename std::vector<int>::const_iterator it_t;

  // initialize+wrong communications
  A.apply_everywhere(c);
  // initialize+communicate
  A.apply_everywhere(c);

  for( int rank = 0; rank < numprocs; ++rank) {
    // test send buffer    
    if( A.get_communicator().rank() == rank) {
      for( it_t it = A.get_communicator().recv_buff()[3].second.begin();
	   it != A.get_communicator().recv_buff()[3].second.end(); ++it)
	ASSERT_EQ(*it, (A.get_communicator().nb()[3].second+1));
    }
    if( A.get_communicator().rank() == rank) {
      for( it_t it = A.get_communicator().recv_buff()[3].first.begin();
	   it != A.get_communicator().recv_buff()[3].first.end(); ++it)
	ASSERT_EQ(*it, (A.get_communicator().nb()[3].first+1));
    }

    MPI_Barrier(MPI_COMM_WORLD);
  }
}

TEST(Comunicator, Unbuffer){
  typedef fields::LocalField<int, 4> intField;
  typedef typename intField::raw_pt raw_pt_t;
  typedef typename std::vector<int> vec_t;
  intField::extents_t e;
  intField::neighbors_t nn;
  std::fill(e.begin(), e.end(), 4);
  intField A(e,1,0,nn);
  cfld<intField> c(A.get_communicator().rank()+1);
  A.apply_everywhere(c);
  A.apply_everywhere(c);
  
  typedef geometry::ZeroBndIterator ZIterator;
  typedef geometry::TBndIterator TIterator;
  

  ZIterator i(A.mk_point(raw_pt_t{0,0,0,0}), e);
  do {
    ASSERT_EQ(A[*i],(A.get_communicator().nb()[3].first+1));
  } while ( (++i).is_good() );

  TIterator f(A.mk_point(raw_pt_t{0,0,0,3}), e);
  do {
    ASSERT_EQ(A[*f],(A.get_communicator().nb()[3].second+1));
  } while ( (++f).is_good() );

}




TEST(Comunicator, Finalize){
  comm::Communicator<fields::LocalField<int, 4> >::finalize();
}

int main(int argc, char **argv) {
  comm::Communicator<fields::LocalField<int, 4> >::init(argc,argv);
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
