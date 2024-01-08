# include(GoogleTest)

# add_executable(test_Geometry test_Geometry.cpp)
# add_executable(test_Point test_Point.cpp)
# target_link_libraries(test_Geometry gtest_main)
# target_link_libraries(test_Point gtest_main)

# gtest_discover_tests(test_Geometry)
# gtest_discover_tests(test_Point)

include(add_unit_test)

add_unit_test(test_Geometry test_Geometry.cpp "" "$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>")
add_unit_test(test_Point test_Point.cpp "" "$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>")