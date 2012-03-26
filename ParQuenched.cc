#include <LocalGluonField.hpp>
#include <Background.h>
#include <Geometry.hpp>
#include <algorithm>
#include <newQCDpt.h>
#include <mpi.h>
#include <iostream>

const int DIM = 4;
const int ORD = 6;

typedef fields::LocalGluonField<bgf::AbelianBgf, ORD, DIM> localGF_t;
typedef localGF_t::neighbors_t nt;
typedef BGptSU3<bgf::AbelianBgf, ORD> ptSU3;

// helper function to conveniently initialize the
// array of neighbors in four dimensions

inline nt neighbors(const int& tu, 
                    const int& td,
                    const int& xu,
                    const int& xd,
                    const int& yu,
                    const int& yd,
                    const int& zu,
                    const int& zd){
  nt result;
  result[0].first = tu;
  result[0].second = td;
  result[1].first = xu;
  result[1].second = xd;
  result[2].first = yu;
  result[2].second = yd;
  result[3].first = zu;
  result[3].second = zd;
  return result;
}


int main(int argc, char *argv[]) {
  // generate an array for to store the lattice extents
  geometry::Geometry<DIM>::extents_t e;
  // we want a L = 4 lattice
  std::fill(e.begin(), e.end(), 4);
  // initialize parallel enviornment
  int rank , numprocs;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (numprocs != 2){
    std::cout << "Please start the test program with two processes!\n";
    return 1;
  }
  // tell the process which are the neighbors
  int zup = (rank == numprocs - 1) ? 0 : rank + 1;
  int zdown = (!rank) ? numprocs - 1 : rank - 1;
  nt n = neighbors(0,0,0,0,0,0, zup, zdown);
  fields::LocalGluonField<bgf::AbelianBgf, 6, DIM> U(e, 1, rank, n);  
  if (rank == 0){
    U.randomize();
    U.test_send_fwd_z();
    // make an iterator for the z-direction, keep z = 4 fixed
    geometry::SliceIterator<DIM> iter = 
      U.mk_slice_iterator(pt::Direction<DIM>(3), 4);
    std::vector<double> traces(ORD + 1, 0);
    while(iter.is_good()){
      pt::Point<DIM> n = iter.yield();
      for (pt::Direction<DIM> mu; mu.is_good(); ++mu){
        traces[0] = U(n,mu).bgf().Tr().re;
        for (int i = 0; i < ORD; ++i) 
          traces[i+1] = U(n,mu)[i].Tr().re;
      }
    }
    std::cout << "real parts traces sent:\n";
    for (int i = 0; i <= ORD; ++i)
      std::cout << i << "\t" << traces[i] << std::endl;
  }
  else {
    U.test_rec_bkw_z();
    // make an iterator for the z-direction, keep z = 0 fixed
    geometry::SliceIterator<DIM> iter = 
      U.mk_slice_iterator(pt::Direction<DIM>(3), 0);
    std::vector<double> traces(ORD + 1, 0);
    while(iter.is_good()){
      pt::Point<DIM> n = iter.yield();
      for (pt::Direction<DIM> mu; mu.is_good(); ++mu){
        traces[0] = U(n,mu).bgf().Tr().re;
        for (int i = 0; i < ORD; ++i) 
          traces[i+1] = U(n,mu)[i].Tr().re;
      }
    }
    std::cout << "real parts of traces received:\n";
    for (int i = 0; i <= ORD; ++i)
      std::cout << i << "\t" << traces[i] << std::endl;
  }
  
  MPI_Finalize();
  return 0;
}
