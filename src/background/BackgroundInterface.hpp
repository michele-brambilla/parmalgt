#ifndef _BACKGROUND_INTERFACE_HPP_
#define _BACKGROUND_INTERFACE_HPP_

#include "BackgroundHelper.hpp"

#include <types/Types.hpp>

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
/// Background field classes.
///
/// Classes to implement trivial and abelian background fields.

namespace bgf {

  
  /// Base class to define interface.
  class BgfBase {
  public:
    /// Left multiplication.
    virtual SU3 ApplyFromLeft ( const SU3 & ) const = 0;
    virtual CVector ApplyFromLeft ( const CVector& ) const = 0;
    /// Right multiplication.
    virtual SU3 ApplyFromRight ( const SU3 & ) const = 0;
    virtual CVector ApplyFromRight ( const CVector & ) const = 0;
    // Add to a SU(3) matrix
    virtual SU3 Add ( const SU3 & ) const = 0;
    // set to zero, or one
    virtual void set_to_one () = 0;
    virtual void set_to_zero () = 0;
    // Trace
    virtual Cplx Tr() const = 0;
    virtual double Norm() const  = 0;
    // Make traceless
    virtual void Trless() = 0;
    // Make (anti-)hermitian and traceless
    virtual void reH() = 0;
  };
};

#endif // _BACKGROUND_INTERFACE_HPP_
