# CMake generated Testfile for 
# Source directory: /mnt/c/Users/Rayan/Research/libami
# Build directory: /mnt/c/Users/Rayan/Research/libami/build
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test([=[fb]=] "test/fb")
set_tests_properties([=[fb]=] PROPERTIES  _BACKTRACE_TRIPLES "/mnt/c/Users/Rayan/Research/libami/CMakeLists.txt;90;add_test;/mnt/c/Users/Rayan/Research/libami/CMakeLists.txt;0;")
add_test([=[helper_functions]=] "test/helper_functions")
set_tests_properties([=[helper_functions]=] PROPERTIES  _BACKTRACE_TRIPLES "/mnt/c/Users/Rayan/Research/libami/CMakeLists.txt;91;add_test;/mnt/c/Users/Rayan/Research/libami/CMakeLists.txt;0;")
add_test([=[num_test]=] "test/num_test")
set_tests_properties([=[num_test]=] PROPERTIES  _BACKTRACE_TRIPLES "/mnt/c/Users/Rayan/Research/libami/CMakeLists.txt;92;add_test;/mnt/c/Users/Rayan/Research/libami/CMakeLists.txt;0;")
add_test([=[construct_tests]=] "test/construct_tests")
set_tests_properties([=[construct_tests]=] PROPERTIES  _BACKTRACE_TRIPLES "/mnt/c/Users/Rayan/Research/libami/CMakeLists.txt;93;add_test;/mnt/c/Users/Rayan/Research/libami/CMakeLists.txt;0;")
subdirs("test")
