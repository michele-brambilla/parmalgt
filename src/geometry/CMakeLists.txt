add_library(Geometry Geometry.cpp Point.cpp)

target_include_directories(Geometry PUBLIC "${PROJECT_BINARY_DIR}"
                                           "${PROJECT_SOURCE_DIR}/src")

add_subdirectory(tests)
