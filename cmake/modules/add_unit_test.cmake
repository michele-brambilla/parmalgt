include(GoogleTest)


# function(add_unit_test test_name test_src_file link_libraries)
#     _add_unit_test(${test_name} ${test_src_file} "NotRoot" "${link_libraries}")
# endfunction()

# function(add_unit_test test_name test_src_file link_libraries)

#     set(test_bin_file test_${test_name})
#     add_executable(${test_bin_file} ${test_src_file})
#     target_link_libraries(${test_bin_file} ${link_libraries} ${CPPUTEST_LIB})
#     add_dependencies(${test_bin_file} cpputest_library)
#     if(${test_type} STREQUAL "AsRoot")
#         target_compile_definitions(${test_bin_file} PUBLIC ASROOT)
#     endif()

#     set(test_commands "cd ${CMAKE_CURRENT_BINARY_DIR} && ${CMAKE_CURRENT_BINARY_DIR}/${test_bin_file}")

#     if(${ENABLE_TESTS_WITH_VALGRIND})
#         set(test_commands "${test_commands} && valgrind --tool=memcheck --leak-check=full --suppressions=valgrind.suppressions --error-exitcode=1 ${CMAKE_CURRENT_BINARY_DIR}/${test_bin_file}")
#     endif()

#     add_test(${test_name} bash -c "${test_commands}")
# endfunction()


function(add_unit_test test_name test_src_file link_libraries include_directories)
    
    add_executable(${test_name} ${test_src_file})

    target_link_libraries(${test_name} gtest_main ${link_libraries})
    target_include_directories(${test_name} PRIVATE ${include_directories})

    gtest_discover_tests(${test_name})

endfunction()