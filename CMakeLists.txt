cmake_minimum_required(VERSION 2.8.12.2)

include(ProcessorCount)
ProcessorCount(N)
if(NOT N EQUAL 0)
    message("Using ${N} cores")
    set(CTEST_BUILD_FLAGS -j${N})
    set(ctest_test_args ${ctest_test_args} PARALLEL_LEVEL ${N})
endif()

project(raptor CXX C)

include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)

set(PROPERTIES CXX_STANDARD 11)
set(PROPERTIES CXX_STANDARD_REQUIRED ON)

add_executable(raptor src/Edge.cpp src/EventQueue.cpp src/GUI.cpp src/Main.cpp src/MessageLogger.cpp src/OctaveWrapper.cpp src/ResourceManager.cpp src/Router.cpp src/Thread.cpp src/Workstation.cpp)

target_include_directories(raptor PUBLIC kshortestpath/include)
target_include_directories(raptor PUBLIC include)

add_subdirectory(kshortestpath)
target_link_libraries(raptor kshortestpath)

if(WIN32)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /EHsc")

	add_definitions(-D_USE_MATH_DEFINES)
endif(WIN32)

if(WIN32)
	include_directories("C:/Users/Tim Hahn/boost_1_66_0")
endif(WIN32)

if(WIN32)
    set(CMAKE_FIND_LIBRARY_PREFIXES "")
    set(CMAKE_FIND_LIBRARY_SUFFIXES ".dll")

	find_library( PTHREAD_VC2_LIB
                  NAMES pthreadVC2.dll
                  HINTS "C:/pthreads/Pre-built.2/dll/x86"
                )

	include_directories("C:/pthreads/Pre-built.2/include")

	target_link_libraries(raptor ${PTHREAD_VC2_LIB})
endif(WIN32)

if(UNIX)
	find_package(Threads REQUIRED)
	target_link_libraries(raptor ${CMAKE_THREAD_LIBS_INIT})
endif(UNIX)

if(WIN32)
    set(CMAKE_FIND_LIBRARY_PREFIXES "")
    set(CMAKE_FIND_LIBRARY_SUFFIXES ".lib" ".dll")

    find_library( OCTAVE_LIB
                  NAMES liboctave-7
                  HINTS "C:/Octave/Octave-5.1.0.0/bin"
                )

    find_library( OCTINTERP_LIB
                  NAMES liboctinterp-7
                  HINTS "C:/Octave/Octave-5.1.0.0/bin"
                )

    include_directories("C:/Octave/Octave-5.1.0.0/mingw32/include/octave-5.1.0/octave")
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
