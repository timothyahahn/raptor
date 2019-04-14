cmake_minimum_required (VERSION 3.14)

project (kshortestpath)

add_library(kshortestpath src/MainP.cpp src/QYDirectedGraph.cpp src/QYKShortestPaths.cpp src/QYShortestPath.cpp)

target_include_directories(kshortestpath PUBLIC include)

