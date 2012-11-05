/// Implementation of the dirac algebra
#include <boost/numeric/ublas/matrix.hpp>
#include "MyMath.h"

namespace Dirac {
  /// Contains the basis of the Dirac algebra
  struct DiracAlgebra {
    typedef boost::numeric::ublas::matrix<Cplx> matrix_t;
    std::vector<matrix_t> basis;
    DiracAlgebra() : basis(6, matrix_t(4,4)) {
      for (int k = 0; k < 5; ++k)
        basis[k].clear();
      Cplx I(0,1);
      // 0) gamma_0
      basis[0](0,2) = 1;
      basis[0](1,3) = 1;
      basis[0](2,0) = 1;
      basis[0](3,1) = 1;
      // 1) gamma_1
      basis[1](0,3) = -I;
      basis[1](1,2) = -I;
      basis[1](2,1) = I;
      basis[1](3,0) = I;
      // 2) gamma_2
      basis[2](0,3) = -1;
      basis[2](1,2) = 1;
      basis[2](2,1) = 1;
      basis[2](3,0) = -1;
      // 3) gamma_3
      basis[3](0,2) = -I;
      basis[3](1,3) = I;
      basis[3](2,0) = I;
      basis[3](3,1) = -I;
      // 4) unit_matrix
      for (int i = 0; i < 4; ++i)
      basis[4](i,i) = 1;
      // 5) gamma_5
      basis[5](0,0) = 1;
      basis[5](1,1) = 1;
      basis[5](2,2) = -1;
      basis[5](3,3) = -1;

    }
    const matrix_t& operator[](const int& i) const {
      return basis[i];
    }
    matrix_t sigma(const int& i, const int& j) const{
      static Cplx Io2(0,0.5);
      return Io2 * boost::numeric::ublas::prod
        (basis[i],basis[j]) - 
        boost::numeric::ublas::prod
        (basis[j],basis[i]);
    }
    matrix_t a_comm(const int& i, const int& j) const{
      return boost::numeric::ublas::prod
        (basis[i],basis[j]) + 
        boost::numeric::ublas::prod
        (basis[j],basis[i]);
    }
  };

  DiracAlgebra::matrix_t dag(const DiracAlgebra::matrix_t &A){
    DiracAlgebra::matrix_t B(A.size1(), A.size2());
    for (int i = 0; i < A.size1(); ++i)
      for (int j = 0; j < A.size2(); ++j){
        B(i,j).re = A(j,i).re;
        B(i,j).im = -A(j,i).im;
      }
    return B;
  }

  DiracAlgebra::matrix_t operator*(const DiracAlgebra::matrix_t &A,
                                   const DiracAlgebra::matrix_t& B){
    return boost::numeric::ublas::prod(A,B);
  }
}
