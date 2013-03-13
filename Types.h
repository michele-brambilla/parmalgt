#ifndef _TYPES_H
#define _TYPES_H
#include <config.h>
#include <complex>
#include <SUN.hpp>
#include <Vector.hpp>

typedef sun::SU<3> SU3;
typedef sun::SU<3>::data_t Cplx;
typedef sun::Vec<3> CVector;

#ifdef HAVE_CXX0X
#include <array>
typedef std::array<Cplx,3> three_vec_t;
template <class C, int n> struct array_t 
{
  typedef std::array<C,n> Type;
};
#elif HAVE_TR1
#include <tr1/array>
typedef std::tr1::array<Cplx,3> three_vec_t;
template <class C, int n> struct array_t 
{
  typedef std::tr1::array<C,n> Type;
};
#endif

#endif
