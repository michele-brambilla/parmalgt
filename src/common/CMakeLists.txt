add_library(Common INTERFACE)

target_include_directories(
  Common
  INTERFACE "$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>")

# add_subdirectory(tests)
