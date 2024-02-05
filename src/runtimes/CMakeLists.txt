add_executable(LocalQuenchSF LocalQuenchSF.cpp)

target_include_directories(
  LocalQuenchSF
  PRIVATE "$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>;$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/deps/sun>;$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/deps/ranlxd>")

target_link_libraries(LocalQuenchSF ranlxd Background Kernels)

add_executable(LocalQuench LocalQuench.cpp)

target_include_directories(
  LocalQuench
  PRIVATE "$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>;$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/deps/sun>;$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/deps/ranlxd>")

target_link_libraries(LocalQuench ranlxd Background Kernels)