include(add_unit_test)

add_unit_test(test_LocalField test_LocalField.cpp "" "$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>;$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/deps/sun>")