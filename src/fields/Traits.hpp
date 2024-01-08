#include <include/cplx.hpp>

#include <complex>
#include <type_traits>

template <typename T> struct undelying_type {
  using type = typename undelying_type<typename T::value_type>::type;
};

template <> struct undelying_type<int> {
  using type = int;
};

template <> struct undelying_type<float> {
  using type = float;
};

template <> struct undelying_type<double> {
  using type = double;
};

template <> struct undelying_type<std::complex<float>> {
  using type = std::complex<float>;
};

template <> struct undelying_type<complex> {
  using type = std::complex<float>;
};