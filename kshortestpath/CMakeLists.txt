cmake_minimum_required (VERSION 2.8.12.2)

project (kshortestpath CXX C)

add_library(kshortestpath src/MainP.cpp src/DirectedGraph.cpp src/KShortestPaths.cpp src/ShortestPath.cpp)

target_include_directories(kshortestpath PUBLIC include)

if(WIN32)
	include_directories("C:/boost_1_70_0")
endif(WIN32)

