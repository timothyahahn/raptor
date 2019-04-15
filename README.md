# raptor

## Linux Installation Instructions

These instructions have been tested on a fresh installation of Fedora 29. Your mileage my vary.

1. Install the required dependancies
```
sudo yum install cmake
sudo yum install gcc-c++
sudo yum install octave-devel
sudo yum install boost-devel
```
2. Run cmake
```
[thahn@localhost raptor]$ cmake .
-- The C compiler identification is GNU 8.3.1
-- The CXX compiler identification is GNU 8.3.1
-- Check for working C compiler: /usr/bin/cc
-- Check for working C compiler: /usr/bin/cc -- works
-- Detecting C compiler ABI info
-- Detecting C compiler ABI info - done
-- Detecting C compile features
-- Detecting C compile features - done
-- Check for working CXX compiler: /usr/bin/c++
-- Check for working CXX compiler: /usr/bin/c++ -- works
-- Detecting CXX compiler ABI info
-- Detecting CXX compiler ABI info - done
-- Detecting CXX compile features
-- Detecting CXX compile features - done
-- Looking for pthread.h
-- Looking for pthread.h - found
-- Looking for pthread_create
-- Looking for pthread_create - not found
-- Check if compiler accepts -pthread
-- Check if compiler accepts -pthread - yes
-- Found Threads: TRUE  
-- Configuring done
-- Generating done
-- Build files have been written to: /home/thahn/raptor
```
3. Build it
```
[thahn@localhost raptor]$ make
Scanning dependencies of target kshortestpath
[  6%] Building CXX object kshortestpath/CMakeFiles/kshortestpath.dir/src/MainP.cpp.o
[ 12%] Building CXX object kshortestpath/CMakeFiles/kshortestpath.dir/src/QYDirectedGraph.cpp.o
[ 18%] Building CXX object kshortestpath/CMakeFiles/kshortestpath.dir/src/QYKShortestPaths.cpp.o
[ 25%] Building CXX object kshortestpath/CMakeFiles/kshortestpath.dir/src/QYShortestPath.cpp.o
[ 31%] Linking CXX static library libkshortestpath.a
[ 31%] Built target kshortestpath
Scanning dependencies of target raptor
[ 37%] Building CXX object CMakeFiles/raptor.dir/src/Edge.cpp.o
[ 43%] Building CXX object CMakeFiles/raptor.dir/src/EventQueue.cpp.o
[ 50%] Building CXX object CMakeFiles/raptor.dir/src/GUI.cpp.o
[ 56%] Building CXX object CMakeFiles/raptor.dir/src/Main.cpp.o
[ 62%] Building CXX object CMakeFiles/raptor.dir/src/MessageLogger.cpp.o
[ 68%] Building CXX object CMakeFiles/raptor.dir/src/OctaveWrapper.cpp.o
[ 75%] Building CXX object CMakeFiles/raptor.dir/src/ResourceManager.cpp.o
[ 81%] Building CXX object CMakeFiles/raptor.dir/src/Router.cpp.o
[ 87%] Building CXX object CMakeFiles/raptor.dir/src/Thread.cpp.o
[ 93%] Building CXX object CMakeFiles/raptor.dir/src/Workstation.cpp.o
[100%] Linking CXX executable raptor
[100%] Built target raptor
```
4. Run it!
