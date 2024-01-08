include(add_unit_test)

add_unit_test(test_PtTypes test_PtTypes.cpp "ranlxd" "$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>;$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/deps/sun>;$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/deps/ranlxd>")