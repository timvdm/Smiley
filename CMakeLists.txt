cmake_minimum_required(VERSION 2.6)
project(Smiley)

# install target for smiley header file
install(FILES smiley.h DESTINATION include)

# examples
add_executable(smiley_print print.cpp)
add_executable(smiley_parse_file parse_file.cpp)

# unit test
enable_testing()
set(TEST_PATH ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
add_executable(smiley_test test.cpp)
add_test(SmileyTest ${TEST_PATH}/smiley_test)
set_tests_properties(SmileyTest PROPERTIES FAIL_REGULAR_EXPRESSION "FAIL")
