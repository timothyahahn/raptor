cmake_minimum_required(VERSION 2.8.12.2)

if(WIN32)
    set(CMAKE_CXX_COMPILER "C:/Program Files/mingw-w64/x86_64-8.1.0-posix-seh-rt_v6-rev0/mingw64/bin/c++.exe")
endif(WIN32)

project(raptor CXX C)

include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)

if(COMPILER_SUPPORTS_CXX11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
else()
        message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()

add_executable(raptor src/Edge.cpp src/EventQueue.cpp src/GUI.cpp src/Main.cpp src/MessageLogger.cpp src/OctaveWrapper.cpp src/ResourceManager.cpp src/Router.cpp src/Thread.cpp src/Workstation.cpp)

target_include_directories(raptor PUBLIC kshortestpath/include)
target_include_directories(raptor PUBLIC include)

add_subdirectory(kshortestpath)
target_link_libraries(raptor kshortestpath)

find_package(Threads REQUIRED)
target_link_libraries(raptor ${CMAKE_THREAD_LIBS_INIT})

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
                )

    find_library( OCTINTERP_LIB
                  NAMES octinterp
                )

	file(GLOB CANDIDATE_OCTAVE_HEADER_DIRS /usr/include/octave*/octave /usr/include/octave)

	find_path( OCTAVE_HEADERS
			   NAMES octave.h
			   HINTS ${CANDIDATE_OCTAVE_HEADER_DIRS}
			 )

    include_directories(${OCTAVE_HEADERS})
endif(UNIX)

if(OCTAVE_LIB)
    target_link_libraries(raptor ${OCTAVE_LIB})
elseif()
    message("OCTAVE_LIB not found, defining NO_OCTAVE")
    add_definitions(-DNO_OCTAVE)
endif(OCTAVE_LIB)

if(OCTINTERP_LIB)
    target_link_libraries(raptor ${OCTINTERP_LIB})
elseif()
    message("OCTINTERP_LIB not found, defining NO_OCTAVE")
    add_definitions(-DNO_OCTAVE)
endif(OCTINTERP_LIB)
