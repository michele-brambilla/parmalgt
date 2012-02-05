#include <iterator>

namespace pta {

  namespace detail {

    template <class T> struct increase_by {
      typename T::const_iterator i;
      typedef typename std::iterator_traits<typename T::iterator>::value_type vt;
      increase_by(const T &obj) : i(obj.begin()) { }
      void operator()(vt& x){ x += *i++; }
    };
  
    template <class T> struct decrease_by {
      typename T::const_iterator i;
      typedef typename std::iterator_traits<typename T::iterator>::value_type vt;
      decrease_by(const T &obj) : i(obj.begin()) { }
      void operator()(vt& x){ x -= *i++; }
    };
  
    template <class T1> struct multiply_by {
      T1 alpha;
      multiply_by (const T1& alph) : alpha(alph) {} 
      template <class T2> void operator()(T2& x) const { x *= alpha; }
    };
  
    template <class T1> struct divide_by {
      T1 alpha;
      divide_by (const T1& alph) : alpha(alph) {} 
      template <class T2> void operator()(T2& x) const { x /= alpha; }
    };

  } // namespace detail

  template <class T> inline detail::increase_by<T> incr (const T& x) { 
    return detail::increase_by<T>(x); 
  }

  template <class T> inline detail::decrease_by<T> decr (const T& x) { 
    return detail::decrease_by<T>(x); 
  }

  template <class T> inline detail::multiply_by<T> mul (const T& x) { 
    return detail::multiply_by<T>(x); 
  }

  template <class T> inline detail::divide_by<T> div (const T& x) { 
    return detail::divide_by<T>(x); 
  }

}
