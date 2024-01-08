add_library(Kernels INTERFACE)

target_include_directories(
  Kernels
  INTERFACE "$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>")

add_subdirectory(tests)
