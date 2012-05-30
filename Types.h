#ifndef _TYPES_H
#define _TYPES_H

#include <config.h>
#include <MyMath.h>

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
