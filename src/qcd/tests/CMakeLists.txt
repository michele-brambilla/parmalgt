include(add_unit_test)

add_unit_test(test_ptSU3 test_ptSU3.cpp "Background;ranlxd" "$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>;$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/deps/sun>;$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/deps/ranlxd>")
add_unit_test(test_BGptSU3 test_BGptSU3.cpp "Background;ranlxd" "$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>;$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/deps/sun>;$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/deps/ranlxd>")