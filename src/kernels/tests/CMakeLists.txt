# include(GoogleTest)

# add_executable(test_Kernels test_Kernels.cpp)

# target_link_libraries(test_Kernels gtest_main)

# target_include_directories(
#   test_Kernels
#   PRIVATE "$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src/kernels>"
#           "$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src/geometry>"
#           "$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src/fields>"
#           "$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src/types>"
#           "$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/deps/sun>")

# gtest_discover_tests(test_Kernels)

include(add_unit_test)

add_unit_test(test_Kernels test_Kernels.cpp "" "$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>;$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/deps/sun>;$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/deps/ranlxd>")