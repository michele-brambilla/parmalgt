#include <config.h>

#ifdef HAVE_CXX0X
#include <array>
typedef std::array<Cplx,3> three_vec_t;
#elif HAVE_TR1
#include <tr1/array>
typedef std::tr1::array<Cplx,3> three_vec_t;
#endif
