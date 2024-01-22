#ifndef UTIL_H
#define UTIL_H

#include <iostream>

#include <string>

// #ifdef _OPENMP
// #include <omp.h>
// #else
namespace kernels {
  inline int omp_get_max_threads() { return 1; }
  inline int omp_get_thread_num() { return 0; }
}
// #endif

#define FLD_INFO(F) \
  typedef typename std_types<F>::ptGluon_t ptGluon;	\
  typedef typename std_types<F>::ptSU3_t ptSU3;		\
  typedef typename std_types<F>::ptsu3_t ptsu3;		\
  typedef typename std_types<F>::bgf_t BGF;		\
  typedef typename std_types<F>::point_t Point;		\
  typedef typename std_types<F>::direction_t Direction;	\
  static const int ORD = std_types<F>::order;		\
  static const int DIM = std_types<F>::n_dim;

namespace util {
  ////////////////////////////////////////////////////////////
  // formated cout for the timings/parameters
  template <typename T>
  inline void pretty_print(const std::string& s, const T& d,
                           const std::string& unit = "",
			   std::ostream& os = std::cout){
    os.width(25); 
    os << s; 
    os.width(0);
    os << ": " << d << unit << std::endl;
  };
}

#endif
