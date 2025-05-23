# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.22

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /mnt/c/Users/Rayan/Research/libami

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /mnt/c/Users/Rayan/Research/libami/build

# Include any dependencies generated for this target.
include test/CMakeFiles/fb.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include test/CMakeFiles/fb.dir/compiler_depend.make

# Include the progress variables for this target.
include test/CMakeFiles/fb.dir/progress.make

# Include the compile flags for this target's objects.
include test/CMakeFiles/fb.dir/flags.make

test/CMakeFiles/fb.dir/fermi_bose.cpp.o: test/CMakeFiles/fb.dir/flags.make
test/CMakeFiles/fb.dir/fermi_bose.cpp.o: ../test/fermi_bose.cpp
test/CMakeFiles/fb.dir/fermi_bose.cpp.o: test/CMakeFiles/fb.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/mnt/c/Users/Rayan/Research/libami/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object test/CMakeFiles/fb.dir/fermi_bose.cpp.o"
	cd /mnt/c/Users/Rayan/Research/libami/build/test && /usr/bin/mpicxx $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT test/CMakeFiles/fb.dir/fermi_bose.cpp.o -MF CMakeFiles/fb.dir/fermi_bose.cpp.o.d -o CMakeFiles/fb.dir/fermi_bose.cpp.o -c /mnt/c/Users/Rayan/Research/libami/test/fermi_bose.cpp

test/CMakeFiles/fb.dir/fermi_bose.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/fb.dir/fermi_bose.cpp.i"
	cd /mnt/c/Users/Rayan/Research/libami/build/test && /usr/bin/mpicxx $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /mnt/c/Users/Rayan/Research/libami/test/fermi_bose.cpp > CMakeFiles/fb.dir/fermi_bose.cpp.i

test/CMakeFiles/fb.dir/fermi_bose.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/fb.dir/fermi_bose.cpp.s"
	cd /mnt/c/Users/Rayan/Research/libami/build/test && /usr/bin/mpicxx $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /mnt/c/Users/Rayan/Research/libami/test/fermi_bose.cpp -o CMakeFiles/fb.dir/fermi_bose.cpp.s

# Object files for target fb
fb_OBJECTS = \
"CMakeFiles/fb.dir/fermi_bose.cpp.o"

# External object files for target fb
fb_EXTERNAL_OBJECTS =

test/fb: test/CMakeFiles/fb.dir/fermi_bose.cpp.o
test/fb: test/CMakeFiles/fb.dir/build.make
test/fb: libami.so
test/fb: lib/libgtest.a
test/fb: lib/libgtest_main.a
test/fb: lib/libgtest.a
test/fb: test/CMakeFiles/fb.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/mnt/c/Users/Rayan/Research/libami/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable fb"
	cd /mnt/c/Users/Rayan/Research/libami/build/test && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/fb.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
test/CMakeFiles/fb.dir/build: test/fb
.PHONY : test/CMakeFiles/fb.dir/build

test/CMakeFiles/fb.dir/clean:
	cd /mnt/c/Users/Rayan/Research/libami/build/test && $(CMAKE_COMMAND) -P CMakeFiles/fb.dir/cmake_clean.cmake
.PHONY : test/CMakeFiles/fb.dir/clean

test/CMakeFiles/fb.dir/depend:
	cd /mnt/c/Users/Rayan/Research/libami/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /mnt/c/Users/Rayan/Research/libami /mnt/c/Users/Rayan/Research/libami/test /mnt/c/Users/Rayan/Research/libami/build /mnt/c/Users/Rayan/Research/libami/build/test /mnt/c/Users/Rayan/Research/libami/build/test/CMakeFiles/fb.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : test/CMakeFiles/fb.dir/depend

