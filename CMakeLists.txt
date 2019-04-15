cmake_minimum_required(VERSION 3.14)

if(WIN32)
	set(CMAKE_C_COMPILER "C:/Program Files/mingw-w64/x86_64-8.1.0-posix-seh-rt_v6-rev0/mingw64/bin/gcc.exe")
	set(CMAKE_CXX_COMPILER "C:/Program Files/mingw-w64/x86_64-8.1.0-posix-seh-rt_v6-rev0/mingw64/bin/g++.exe")
endif(WIN32)

project(raptor)

set(CMAKE_SUPPRESS_REGENERATION true)

add_executable(raptor src/Edge.cpp src/EventQueue.cpp src/GUI.cpp src/Main.cpp src/MessageLogger.cpp src/OctaveWrapper.cpp src/ResourceManager.cpp src/Router.cpp src/Thread.cpp src/Workstation.cpp)

target_include_directories(raptor PUBLIC kshortestpath/include)
target_include_directories(raptor PUBLIC include)

add_subdirectory(kshortestpath)
target_link_libraries(raptor kshortestpath)

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)
target_link_libraries(raptor Threads::Threads)

if(WIN32)
	include_directories("C:/Users/Tim Hahn/boost_1_66_0")
endif(WIN32)

if(WIN32)
	SET(CMAKE_FIND_LIBRARY_PREFIXES "")
	SET(CMAKE_FIND_LIBRARY_SUFFIXES ".lib" ".dll")

    find_library( OCTAVE_LIB
                  NAMES liboctave-4
                  HINTS "C:/Octave/Octave-4.2.2/bin"
                )

    find_library( OCTINTERP_LIB
                  NAMES liboctinterp-4
                  HINTS "C:/Octave/Octave-4.2.2/bin"
                )
                
	include_directories("C:/Octave/Octave-4.2.2/include")
endif(WIN32)

if(UNIX)
	SET(CMAKE_FIND_LIBRARY_PREFIXES "lib")
	SET(CMAKE_FIND_LIBRARY_SUFFIXES ".so" ".a")

    find_library( OCTAVE_LIB
                  NAMES octave
                  HINTS "/usr/lib64/octave/4.2.2"
                )

    find_library( OCTINTERP_LIB
                  NAMES octinterp
                  HINTS "/usr/lib64/octave/4.2.2"
                )

	include_directories(/usr/include/octave-4.2.2/octave)
endif(UNIX)

target_link_libraries(raptor ${OCTAVE_LIB})
target_link_libraries(raptor ${OCTINTERP_LIB})
