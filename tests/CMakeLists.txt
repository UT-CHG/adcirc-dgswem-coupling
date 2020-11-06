cmake_minimum_required(VERSION 3.6)

#------------------------------------------------------------------------------#
# List all files containing tests.
list(APPEND tests_sources
    test_adcoupler_main.F95
)
list(REMOVE_DUPLICATES tests_sources)

list(APPEND tests_sources test_main.F95) # Always include test_main.F95
list(APPEND lib_depends ${LIBRARY_NAME})

set(TEST_MAIN unit_tests)   # Default name for test executable (change if you wish).
set(TEST_RUNNER_PARAMS "")  # Any arguemnts to feed the test runner (change as needed).

#------------------------------------------------------------------------------#
#                         Make Tests (no change needed).
#------------------------------------------------------------------------------#
add_executable(${TEST_MAIN} ${tests_sources})
target_link_libraries(${TEST_MAIN} PRIVATE ${lib_depends})
set_target_properties(${TEST_MAIN} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

set_target_properties(${TEST_MAIN} PROPERTIES 
    Fortran_STANDARD 08 
    Fortran_STANDARD_REQUIRED YES 
    Fortran_EXTENSIONS NO
)

add_test(
    NAME ${LIBRARY_NAME}.${TEST_MAIN}
    COMMAND ${TEST_MAIN} ${TEST_RUNNER_PARAMS}
)

#------------------------------------------------------------------------------#
# The following is inspired from https://github.com/jhbell/cmake-gcov
# Adds a 'coverage' target.
if (CMAKE_BUILD_TYPE STREQUAL "Coverage")
    # Req.d for getting *.o instead of *.F95.o obj files during cmake building
    set(CMAKE_Fortran_OUTPUT_EXTENSION_REPLACE ON)

    # Create the gcov target. Run coverage tests with 'make gcov'
    set(TEST_MAIN_OBJ_DIR "${CMAKE_BINARY_DIR}/tests/CMakeFiles/${TEST_MAIN}.dir")
    add_custom_target(coverage
        COMMAND mkdir -p coverage
        COMMAND ${CMAKE_MAKE_PROGRAM} test
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
        )
    add_custom_command(TARGET coverage
        COMMAND echo "=================== GCOV ===================="
        COMMAND rm -f ${TEST_MAIN_OBJ_DIR}/test_main.gcno  ${TEST_MAIN_OBJ_DIR}/test_main.gcda  ${TEST_MAIN_OBJ_DIR}/test_main.o
        COMMAND gcov -b ${CMAKE_SOURCE_DIR}/tests/*.F95 -o ${TEST_MAIN_OBJ_DIR} > CoverageSummary.txt
        COMMAND cat CoverageSummary.txt
        COMMAND echo "-- Coverage files have been output to ${CMAKE_BINARY_DIR}/coverage"
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/coverage  # Need separate command for this line
        )

    add_dependencies(coverage ${TEST_MAIN})

    # Make sure to clean up the coverage folder when running make clean
    set_property(DIRECTORY PROPERTY ADDITIONAL_MAKE_CLEAN_FILES 
        "${CMAKE_BINARY_DIR}/coverage"
        )
    # Make sure to clean up the gcov files when running make clean
    set_property(DIRECTORY PROPERTY ADDITIONAL_MAKE_CLEAN_FILES 
        "${TEST_MAIN_OBJ_DIR}/main.gcno"
        "${TEST_MAIN_OBJ_DIR}/main.gcda"
        )

endif() #CMAKE_BUILD_TYPE STREQUAL "Coverage"

