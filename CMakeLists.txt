cmake_minimum_required(VERSION 3.14)

project(raptor)

add_subdirectory(kshortestpath)

add_executable(raptor src/Edge.cpp src/EventQueue.cpp src/GUI.cpp src/Main.cpp src/MessageLogger.cpp src/ResourceManager.cpp src/Router.cpp src/Thread.cpp src/Workstation.cpp)

target_include_directories(raptor PUBLIC kshortestpath/include)
target_include_directories(raptor PUBLIC include)
