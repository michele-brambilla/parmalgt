add_library(Kernels)

target_sources(
  Kernels
  PRIVATE IO.cpp
  INTERFACE IO.hpp Kernels.hpp Plaquette.hpp Project.hpp RungeKutta.hpp util.hpp uparam.hpp)

target_include_directories(
  Kernels
  PRIVATE "$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>;$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/deps/sun>")

add_subdirectory(tests)
