add_library(Fields INTERFACE)

target_include_directories(
  Fields INTERFACE "$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>")

add_subdirectory(tests)
