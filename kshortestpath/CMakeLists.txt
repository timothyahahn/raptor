cmake_minimum_required (VERSION 3.14)

project (kshortestpath)

add_library(kshortestpath src/MainP.cpp src/QYDirectedGraph.cpp src/QYKShortestPaths.cpp src/QYShortestPath.cpp)

target_include_directories(kshortestpath PUBLIC include)

if(WIN32)
	include_directories("C:\\Users\\Tim Hahn\\boost_1_66_0")
endif(WIN32)