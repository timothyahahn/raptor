cmake_minimum_required(VERSION 2.8.12.2)

macro(use_cxx11)
    if (CMAKE_VERSION VERSION_LESS "3.1")
        if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --std=c++11")
        endif()
    else()
        set(CMAKE_CXX_STANDARD 11)
    endif()
endmacro(use_cxx11)

include(ProcessorCount)
ProcessorCount(N)
if(NOT N EQUAL 0)
    message("Using ${N} cores")
    set(CTEST_BUILD_FLAGS -j${N})
    set(ctest_test_args ${ctest_test_args} PARALLEL_LEVEL ${N})
endif()

project(raptor CXX C)

set(CMKAE_BUILD_TYPE RelWithDebInfo)

use_cxx11()

add_executable(raptor src/Edge.cpp src/EventQueue.cpp src/GUI.cpp src/Main.cpp src/MessageLogger.cpp src/OctaveWrapper.cpp src/ResourceManager.cpp src/Router.cpp src/Thread.cpp)

target_include_directories(raptor PUBLIC kshortestpath/include)
target_include_directories(raptor PUBLIC include)

add_subdirectory(kshortestpath)
target_link_libraries(raptor kshortestpath)

if(WIN32)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /EHsc")

    add_definitions(-D_USE_MATH_DEFINES)
endif(WIN32)

if(WIN32)
    set(CMAKE_FIND_LIBRARY_PREFIXES "")
    set(CMAKE_FIND_LIBRARY_SUFFIXES ".lib")

	find_library( PTHREAD_VC2_LIB
                  NAMES pthreadVC2.lib
                  HINTS "C:/pthreads/Pre-built.2/lib/x64"
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

    file(GLOB CANDIDATE_OCTAVE_LIB_DIRS "C:/Octave/Octave-5.1.0.0/bin")

    find_library( OCTAVE_LIB
                  NAMES liboctave-7
                  HINTS ${CANDIDATE_OCTAVE_LIB_DIRS}
                )

    find_library( OCTINTERP_LIB
                  NAMES liboctinterp-7
                  HINTS ${CANDIDATE_OCTAVE_LIB_DIRS}
                )

    if(NOT OCTAVE_LIB OR NOT OCTINTERP_LIB)
        message("Warning: Octave libraries not found in: ${CANDIDATE_OCTAVE_LIB_DIRS}")
        message("Defining NO_OCTAVE")
        add_definitions(-DNO_OCTAVE)
    else()
        target_link_libraries(raptor ${OCTAVE_LIB})
        target_link_libraries(raptor ${OCTINTERP_LIB})
    endif()

    file(GLOB CANDIDATE_OCTAVE_HEADER_DIRS "C:/Octave/Octave-5.1.0.0/mingw32/include/octave-5.1.0")

    find_path( OCTAVE_HEADERS
               NAMES octave/octave.h
               HINTS ${CANDIDATE_OCTAVE_HEADER_DIRS}
	     )

    if(NOT OCTAVE_HEADERS)
        message("Warning: Octave headers not found in: ${CANDIDATE_OCTAVE_HEADER_DIRS}")
        message("Defining NO_OCTAVE")
        add_definitions(-DNO_OCTAVE)
    else()
        include_directories(${OCTAVE_HEADERS})
    endif()

endif(WIN32)

if(UNIX)
    SET(CMAKE_FIND_LIBRARY_PREFIXES "lib")
    SET(CMAKE_FIND_LIBRARY_SUFFIXES ".so" ".a")

    file(GLOB CANDIDATE_OCTAVE_LIB_DIRS /usr/lib64/octave/* /usr/lib64)

    find_library( OCTAVE_LIB
                  NAMES octave
                  HINTS ${CANDIDATE_OCTAVE_LIB_DIRS}
                )

    find_library( OCTINTERP_LIB
                  NAMES octinterp
                  HINTS ${CANDIDATE_OCTAVE_LIB_DIRS}
                )

    if(NOT OCTAVE_LIB OR NOT OCTINTERP_LIB)
        message("Warning: Octave libraries not found in: ${CANDIDATE_OCTAVE_LIB_DIRS}")
        message("Defining NO_OCTAVE")
        add_definitions(-DNO_OCTAVE)
    else()
        target_link_libraries(raptor ${OCTAVE_LIB})
        target_link_libraries(raptor ${OCTINTERP_LIB})
    endif()

    file(GLOB CANDIDATE_OCTAVE_HEADER_DIRS /usr/include/octave* /usr/include/octave*/octave /usr/include)

    find_path( OCTAVE_HEADERS
               NAMES octave/octave.h
               HINTS ${CANDIDATE_OCTAVE_HEADER_DIRS}
	     )

    find_path( OCTAVE_SUB_HEADERS
	       NAMES str-vec.h
	       HINTS ${CANDIDATE_OCTAVE_HEADER_DIRS}
	     )

     if(NOT OCTAVE_HEADERS OR NOT OCTAVE_SUB_HEADERS)
        message("Warning: Octave headers not found in: ${CANDIDATE_OCTAVE_HEADER_DIRS}")
        message("Defining NO_OCTAVE")
        add_definitions(-DNO_OCTAVE)
    else()
        include_directories(${OCTAVE_HEADERS})
	include_directories(${OCTAVE_SUB_HEADERS})
    endif()
endif(UNIX)

if(WIN32)
    set(CMAKE_FIND_LIBRARY_PREFIXES "")
    set(CMAKE_FIND_LIBRARY_SUFFIXES ".lib" ".dll")

    file(GLOB CANDIDATE_ALLEGRO_LIB_DIRS "C:/Octave/Octave-5.1.0.0/bin")

    find_library( ALLEGRO_LIB
                  NAMES allegro
                  HINTS CANDIDATE_ALLEGRO_LIB_DIRS
                )

    if(NOT ALLEGRO_LIB)
        message("Warning: Allegro libraries not found in: ${CANDIDATE_ALLEGRO_LIB_DIRS}")
        message("Defining NO_ALLEGRO")
        add_definitions(-DNO_ALLEGRO)
    else()
        target_link_libraries(raptor ${ALLEGRO_LIB})
    endif()

    file(GLOB CANDIDATE_ALLEGRO_HEADER_DIRS "C:/allegro5/include")

    find_path( ALLEGRO_HEADERS
               NAMES allegro5/allegro.h
               HINTS ${CANDIDATE_ALLEGRO_HEADER_DIRS}
	     )

    if(NOT ALLEGRO_HEADERS)
        message("Warning: Allegro headers not found in: ${CANDIDATE_ALLEGRO_HEADER_DIRS}")
        message("Defining NO_ALLEGRO")
        add_definitions(-DNO_ALLEGRO)
    else()
        include_directories(${ALLEGRO_HEADERS})
    endif()
endif(WIN32)

if(UNIX)
    SET(CMAKE_FIND_LIBRARY_PREFIXES "lib")
    SET(CMAKE_FIND_LIBRARY_SUFFIXES ".so" ".a")

    file(GLOB CANDIDATE_ALLEGRO_LIB_DIRS /usr/lib64)

    find_library( ALLEGRO_LIB
                  NAMES allegro
                  HINTS ${CANDIDATE_ALLEGRO_LIB_DIRS}
                )

    if(NOT ALLEGRO_LIB)
        message("Warning: Allegro libraries not found in: ${CANDIDATE_ALLEGRO_LIB_DIRS}")
        message("Defining NO_ALLEGRO")
        add_definitions(-DNO_ALLEGRO)
    else()
        target_link_libraries(raptor ${ALLEGRO_LIB})
    endif()

    file(GLOB CANDIDATE_ALLEGRO_HEADER_DIRS /usr/include/allegro5 /usr/include)

    find_path( ALLEGRO_HEADERS
               NAMES allegro5/allegro.h
               HINTS ${CANDIDATE_ALLEGRO_HEADER_DIRS}
	     )

    if(NOT ALLEGRO_HEADERS)
        message("Warning: Allegro headers not found in: ${CANDIDATE_ALLEGRO_HEADER_DIRS}")
        message("Defining NO_ALLEGRO")
        add_definitions(-DNO_ALLEGRO)
    else()
        include_directories(${ALLEGRO_HEADERS})
    endif()
endif(UNIX)

