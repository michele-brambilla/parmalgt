add_library(Background)
target_sources(
  Background
  PRIVATE AbelianBackground.cpp
  INTERFACE BackgroundHelper.hpp BackgroundInterface.hpp AbelianBackground.hpp TrivialBackground.hpp)

target_include_directories(
  Background
  PRIVATE "$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>"
          "$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/deps/sun>"
          "$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/deps/ranlxd>")

add_subdirectory(tests)
