# include(GoogleTest)

# add_executable(test_BackgroundHelper test_BackgroundHelper.cpp)

# target_link_libraries(test_BackgroundHelper gtest_main)

# add_executable(test_TrivialBackground test_TrivialBackground.cpp)

# target_link_libraries(test_TrivialBackground gtest_main)

# # target_include_directories(
# #   test_BackgroundHelper
# #   PRIVATE "$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>"
# #           "$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/deps/sun>")
# target_include_directories(
#   test_TrivialBackground
#   PRIVATE "$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>"
#           "$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/deps/sun>")

# gtest_discover_tests(test_BackgroundHelper)
# gtest_discover_tests(test_TrivialBackground)


include(add_unit_test)

add_unit_test(test_BackgroundHelper test_BackgroundHelper.cpp "" "$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>;$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/deps/sun>")
add_unit_test(test_TrivialBackground test_TrivialBackground.cpp "" "$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>;$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/deps/sun>")