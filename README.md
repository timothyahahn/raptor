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
2. Run cmake (you will need to update the paths based on the versions installed and their location)
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
4. Run it with no parameters
```
[thahn@localhost raptor]$ ./raptor 
Usage: ./raptor <Topology> <Wavelengths> <Random Seed> <Thread Count> <Iteration Count> <Probe Count>

Topology: NSF, Mesh, Mesh6x6, Mesh8x8, Mesh10x10
Wavelengths: 21, 41, 81, 161, 321, 641, 1281
Random Seed: <any valid unsigned short int>
Thread Count: <maximum number of threads to run, 1 to n>
Iteration Count: <number of iterations, 1 to n>
Probe Count: <number of probes, 1 to n>
```
5. Run it for real!
```
[thahn@localhost raptor]$ ./raptor NSF 21 123456789 4 10 10
21:26:26 [] Reading Quality Parameters from input/Quality-NSF-21.txt file.
21:26:26 [] 	arrival_interval = 250.000000
21:26:26 [] 	duration = 250.000000
21:26:26 [] 	nonlinear_halfwin = 10
21:26:26 [] 	halfwavelength = 10
21:26:26 [] 	sys_wavelength = 21
21:26:26 [] 	fc = 193100001574912.000000
21:26:26 [] 	f_step = 49999998976.000000
21:26:26 [] 	channel_power = 0.001000
21:26:26 [] 	L = 80.000000
21:26:26 [] 	alphaDB = 0.250000
21:26:26 [] 	D = 0.004000
21:26:26 [] 	S = 80000.000000
21:26:26 [] 	gamma = 2.600000
21:26:26 [] 	QFactor_factor = 0.950000
21:26:26 [] 	EDFA_Noise_Figure = 3.500000
21:26:26 [] 	EDFA_Gain = 21.600000
21:26:26 [] 	B_w = 1.000000e+10
21:26:26 [] 	usage_update_interval = 4.000000
21:26:26 [] 	gui_update_interval = 6
21:26:26 [] 	beta = 1.200000
21:26:26 [] 	refractive_index = 1.500000
21:26:26 [] 	q_factor_stats = 1
21:26:26 [] 	detailed_log = 0
21:26:26 [] 	DP_alpha = 1.000000
21:26:26 [] 	ACO_ants = 20
21:26:26 [] 	ACO_alpha = 1.000000
21:26:26 [] 	ACO_beta = 5.000000
21:26:26 [] 	ACO_rho = 0.100000
21:26:26 [] 	MM_ACO_gamma = 0.100000
21:26:26 [] 	MM_ACO_N_iter = 20
21:26:26 [] 	MM_ACO_N_reset = 2
Calculating router distances...done.
21:26:26 [] Created 4 threads.

21:26:55 [] **ALGORITHM = DP-Q-FF, WORKS = 280, PROBE = PARALLEL, QA = 1**
21:26:55 [] OVERALL BLOCKING (658/4438) = 0.148265
21:26:55 [] COLLISIONS (0/4438) = 0.000000
21:26:55 [] BAD QUALITY (0/4438) = 0.000000
21:26:55 [] NON RESOURCES (658/4438) = 0.148265
21:26:55 [] AVERAGE PROBES PER REQUEST (4438/4438) = 1.000000
21:26:55 [] AVERAGE REQUEST DELAY TIME (67.459469/3780) = 0.017846
21:26:55 [] AVERAGE CONNECTION HOP COUNT (9209/3780) = 2.436243
21:26:55 [] AVERAGE CONNECTION SPAN COUNT (133974/3780) = 35.442856
21:26:55 [] AVERAGE ASE NOISE (0.000110/3780) = 2.913904e-08
21:26:55 [] AVERAGE FWM NOISE (0.000022/3780) = 5.922165e-09
21:26:55 [] AVERAGE XPM NOISE (0.000000/3780) = 0.000000e+00
21:26:55 [] AVERAGE RA RUN TIME (14.000000/4438) = 3.154574e-03
21:26:55 [] DROPPED CONNECTIONS (597/4438) = 0.134520
21:26:55 [] OVERALL W/DROPPED (1255/4438) = 0.282785
21:26:55 [] INITIAL Q: MIN = 5.860601, MAX = 12.420026, AVG = 7.720318
21:26:55 [] AVERAGE Q: MIN = 3.319328, MAX = 12.417975, AVG = 7.607907
21:26:55 [] % TIME Q BELOW: MIN = 0.000000, MAX = 0.998599, AVG = 0.086621
21:26:55 [] ***********************************************
```
## Windows Installation Instructions

These instructions have been tested on a fresh installation of Windows 10. This is a work in progress...

1. Install the required dependancies.
* boost - https://www.boost.org/users/download/
* cmake - https://cmake.org/download/
* mingw-w64 - https://mingw-w64.org/doku.php/download/mingw-builds
* octave - https://www.gnu.org/software/octave/download.html
2. Run cmake (you will need to update the paths based on the versions installed and their location)
```
C:\Users\Tim Hahn\raptor>cmake . -G "MinGW Makefiles"
-- The C compiler identification is GNU 8.1.0
-- The CXX compiler identification is GNU 8.1.0
-- Check for working C compiler: C:/Program Files/mingw-w64/x86_64-8.1.0-posix-seh-rt_v6-rev0/mingw64/bin/gcc.exe
-- Check for working C compiler: C:/Program Files/mingw-w64/x86_64-8.1.0-posix-seh-rt_v6-rev0/mingw64/bin/gcc.exe -- works
-- Detecting C compiler ABI info
-- Detecting C compiler ABI info - done
-- Detecting C compile features
-- Detecting C compile features - done
-- Check for working CXX compiler: C:/Program Files/mingw-w64/x86_64-8.1.0-posix-seh-rt_v6-rev0/mingw64/bin/g++.exe
-- Check for working CXX compiler: C:/Program Files/mingw-w64/x86_64-8.1.0-posix-seh-rt_v6-rev0/mingw64/bin/g++.exe -- works
-- Detecting CXX compiler ABI info
-- Detecting CXX compiler ABI info - done
-- Detecting CXX compile features
-- Detecting CXX compile features - done
-- Looking for pthread.h
-- Looking for pthread.h - found
-- Looking for pthread_create
-- Looking for pthread_create - found
-- Found Threads: TRUE
-- Configuring done
-- Generating done
-- Build files have been written to: C:/Users/Tim Hahn/raptor
```
3. Build the program
```
C:\Users\Tim Hahn\raptor>mingw32-make
Scanning dependencies of target kshortestpath
[  6%] Building CXX object kshortestpath/CMakeFiles/kshortestpath.dir/src/MainP.cpp.obj
[ 12%] Building CXX object kshortestpath/CMakeFiles/kshortestpath.dir/src/QYDirectedGraph.cpp.obj
[ 18%] Building CXX object kshortestpath/CMakeFiles/kshortestpath.dir/src/QYKShortestPaths.cpp.obj
[ 25%] Building CXX object kshortestpath/CMakeFiles/kshortestpath.dir/src/QYShortestPath.cpp.obj
[ 31%] Linking CXX static library libkshortestpath.a
[ 31%] Built target kshortestpath
Scanning dependencies of target raptor
[ 37%] Building CXX object CMakeFiles/raptor.dir/src/Edge.cpp.obj
In file included from C:/PROGRA~1/MINGW-~1/X86_64~1.0-P/mingw64/lib/gcc/x86_64-w64-mingw32/8.1.0/include/c++/system_error:39,
                 from C:/PROGRA~1/MINGW-~1/X86_64~1.0-P/mingw64/lib/gcc/x86_64-w64-mingw32/8.1.0/include/c++/bits/ios_base.h:46,
                 from C:/PROGRA~1/MINGW-~1/X86_64~1.0-P/mingw64/lib/gcc/x86_64-w64-mingw32/8.1.0/include/c++/ios:42,
                 from C:/PROGRA~1/MINGW-~1/X86_64~1.0-P/mingw64/lib/gcc/x86_64-w64-mingw32/8.1.0/include/c++/ostream:38,
                 from C:/PROGRA~1/MINGW-~1/X86_64~1.0-P/mingw64/lib/gcc/x86_64-w64-mingw32/8.1.0/include/c++/iostream:39,
                 from C:/Users/TIMHAH~1/raptor/include/Thread.h:30,
                 from C:\Users\Tim Hahn\raptor\src\Edge.cpp:23:
C:/PROGRA~1/MINGW-~1/X86_64~1.0-P/mingw64/lib/gcc/x86_64-w64-mingw32/8.1.0/include/c++/x86_64-w64-mingw32/bits/error_constants.h:53:25: error: 'EBADMSG' was not declared in this scope
       bad_message =     EBADMSG,
```