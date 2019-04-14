cmake_minimum_required(VERSION 3.14)

project(raptor)

add_executable(raptor src/Edge.cpp src/EventQueue.cpp src/GUI.cpp src/Main.cpp src/MessageLogger.cpp src/OctaveWrapper.cpp src/ResourceManager.cpp src/Router.cpp src/Thread.cpp src/Workstation.cpp)

target_include_directories(raptor PUBLIC kshortestpath/include)
target_include_directories(raptor PUBLIC include)

add_subdirectory(kshortestpath)
target_link_libraries(raptor kshortestpath)

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)
target_link_libraries(raptor Threads::Threads)

find_library(OCTAVE_LIB octave HINTS /usr/lib64/octave/4.2.2)
target_link_libraries(raptor "${OCTAVE_LIB}")

find_library(OCTINTERP_LIB octinterp HINTS /usr/lib64/octave/4.2.2)
target_link_libraries(raptor "${OCTINTERP_LIB}")

include_directories(/usr/include/octave-4.2.2/octave)
