add_library(Methods INTERFACE)

target_include_directories(
  Methods
  INTERFACE "$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>")

add_subdirectory(tests)
