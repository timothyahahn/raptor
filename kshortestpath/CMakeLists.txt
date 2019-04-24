cmake_minimum_required (VERSION 2.8.12.2)

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

project (kshortestpath CXX)

use_cxx11()

add_library(kshortestpath src/MainP.cpp src/DirectedGraph.cpp src/KShortestPaths.cpp src/ShortestPath.cpp)

target_include_directories(kshortestpath PUBLIC include)

if(WIN32)
	include_directories("C:/boost_1_70_0")
endif(WIN32)

