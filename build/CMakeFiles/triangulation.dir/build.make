# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.16

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
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
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/t225/project_2

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/t225/project_2/build

# Include any dependencies generated for this target.
include CMakeFiles/triangulation.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/triangulation.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/triangulation.dir/flags.make

CMakeFiles/triangulation.dir/src/main.cpp.o: CMakeFiles/triangulation.dir/flags.make
CMakeFiles/triangulation.dir/src/main.cpp.o: ../src/main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/t225/project_2/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/triangulation.dir/src/main.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/triangulation.dir/src/main.cpp.o -c /home/t225/project_2/src/main.cpp

CMakeFiles/triangulation.dir/src/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/triangulation.dir/src/main.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/t225/project_2/src/main.cpp > CMakeFiles/triangulation.dir/src/main.cpp.i

CMakeFiles/triangulation.dir/src/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/triangulation.dir/src/main.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/t225/project_2/src/main.cpp -o CMakeFiles/triangulation.dir/src/main.cpp.s

CMakeFiles/triangulation.dir/src/triangulation.cpp.o: CMakeFiles/triangulation.dir/flags.make
CMakeFiles/triangulation.dir/src/triangulation.cpp.o: ../src/triangulation.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/t225/project_2/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/triangulation.dir/src/triangulation.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/triangulation.dir/src/triangulation.cpp.o -c /home/t225/project_2/src/triangulation.cpp

CMakeFiles/triangulation.dir/src/triangulation.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/triangulation.dir/src/triangulation.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/t225/project_2/src/triangulation.cpp > CMakeFiles/triangulation.dir/src/triangulation.cpp.i

CMakeFiles/triangulation.dir/src/triangulation.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/triangulation.dir/src/triangulation.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/t225/project_2/src/triangulation.cpp -o CMakeFiles/triangulation.dir/src/triangulation.cpp.s

# Object files for target triangulation
triangulation_OBJECTS = \
"CMakeFiles/triangulation.dir/src/main.cpp.o" \
"CMakeFiles/triangulation.dir/src/triangulation.cpp.o"

# External object files for target triangulation
triangulation_EXTERNAL_OBJECTS =

triangulation: CMakeFiles/triangulation.dir/src/main.cpp.o
triangulation: CMakeFiles/triangulation.dir/src/triangulation.cpp.o
triangulation: CMakeFiles/triangulation.dir/build.make
triangulation: libCGAL_Qt5_moc_and_resources.a
triangulation: /usr/lib/x86_64-linux-gnu/libgmpxx.so
triangulation: /usr/lib/x86_64-linux-gnu/libmpfr.so
triangulation: /usr/lib/x86_64-linux-gnu/libgmp.so
triangulation: /usr/lib/x86_64-linux-gnu/libQt5OpenGL.so.5.12.8
triangulation: /usr/lib/x86_64-linux-gnu/libQt5Widgets.so.5.12.8
triangulation: /usr/lib/x86_64-linux-gnu/libQt5Gui.so.5.12.8
triangulation: /usr/lib/x86_64-linux-gnu/libQt5Core.so.5.12.8
triangulation: CMakeFiles/triangulation.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/t225/project_2/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX executable triangulation"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/triangulation.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/triangulation.dir/build: triangulation

.PHONY : CMakeFiles/triangulation.dir/build

CMakeFiles/triangulation.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/triangulation.dir/cmake_clean.cmake
.PHONY : CMakeFiles/triangulation.dir/clean

CMakeFiles/triangulation.dir/depend:
	cd /home/t225/project_2/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/t225/project_2 /home/t225/project_2 /home/t225/project_2/build /home/t225/project_2/build /home/t225/project_2/build/CMakeFiles/triangulation.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/triangulation.dir/depend
