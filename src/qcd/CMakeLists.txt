add_library(qcd INTERFACE)

target_include_directories(
  Types INTERFACE "$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>;$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/deps/sun>;$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/deps/ranlxd>")

add_subdirectory(tests)
