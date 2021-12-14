[![CreamyStrand](http://www.cs.columbia.edu/cg/raymond/pasta_banner.jpg)](http://www.cs.columbia.edu/cg/creamystrand/)

CreamyStrand
================
CreamyStrand is an open source project for the physical simulation of the coupling between hairs and shear-dependent liquid. It supports Linux and Windows, and licensed under the Mozilla Public License v. 2.0.

We would like to hear from you if you appreciate this work.

It is the original implementation of paper *A Multi-Scale Model for Coupling Strands with Shear-Dependent Liquid* ( refer to our project page for more details: http://www.cs.columbia.edu/cg/creamystrand/ ). This code base contains the following parts:

 - A strand simulator adapted from the code of ADONIS (http://www.cs.columbia.edu/cg/adonis/), which adopts discrete elastic rods and nonlinear integration to simulate hairs.
 - A bulk liquid simulator for both shear-dependent and Newtonian liquid, discretized with augmented, moving least sqaures material point method (AMLS-MPM).
 - A reduced-dimensional flow simulator that handles shear-dependent liquid on a strand's surface.
 - A framework coupling the dynamics between strands, the bulk liquid, and reduced-dimensional flows.

Dependencies
--------------------
CreamyStrand depends on following libraries (some of them have been included in the code base, marked with asterisk):

- Eigen* (http://eigen.tuxfamily.org/)
- RapidXML* (http://rapidxml.sourceforge.net)
- tclap* (http://tclap.sourceforge.net)
- libIGL* (https://github.com/libigl/libigl)
- ANN* (https://www.cs.umd.edu/~mount/ANN/)
- So-Bogus* (http://gdaviet.fr/code/bogus/)
- Intel MKL (https://software.intel.com/en-us/parallel-studio-xe/choose-download)
- GLUT (https://www.opengl.org/resources/libraries/glut/)
- Boost 1.58+ (https://www.boost.org/)

On Linux-based systems, most of the dependencies are either included, or can be easily installed with the APT package handling utility. For Intel MKL, you may download and install from the link provided above. For GLUT on Windows, you may acquire the binaries from the NVIDIA Cg Toolkit (https://developer.nvidia.com/cg-toolkit)

For parallelization, CreamyStrand uses the OpenMP library. 

On Windows you may need manually download and compile some of the dependencies. More details refer to the compilation section below.

Compilation
-----------------
CreamyStrand has been tested with GCC 4.8+ (under Linux), and MSVC (under Windows 10 and Visual Studio 2019).

To compile CreamyStrand, you'll need CMake or CMake-GUI (https://cmake.org).

Before the following process, you'd make sure the environment variable `BOOST_ROOT` has been set to the installation location of Boost (the directory that contains `include` and `lib` folders of Boost).

Command Line:
1. make a directory, say, *build*, with *mkdir build*, enter the *build* directory, type *cmake ..*
2. Optionally you can adjust the options with *ccmake ..* In some cases there can be some packages that cmake cannot find. You need to manually specify their paths through ccmake then.
3. type *make* to compile the code. For speeding up the compilation process you may use *make -j*.

CMake-GUI:
1. open CMake-GUI, enter the correct directory for source code and build. Then click *Configure*, choose your installed version of the compiler (on Windows, choose the correct version of Microsoft Visual Studio).
2. after configuration you may find several libraries not found (with notifications of errors), check the *Advanced* box and *specify those missing header path and libraries manually*. For example, if Eigen is missing, then please specify the EIGEN3_INCLUDE_DIR to the path of directory we provided. For the ones we have not provided, you need to download and compile them, and then specify the missing directories to the path containing your headers or compiled libraries. Please make sure you have picked the libraries corresponding to the architecture you have selected (say, 32-bit libraries for x86, and 64-bit libraries for x64).
3. click generate after fixing all missing variables to generate the makefile (on Windows, the generated file is a Visual Studio solution).
4. compile the code.
5. On Windows, before running the demo, all the binary dynamic linking libraries (DLLs) for your dependencies (especially, the Intel MKL/runtime libraries) should be accessible through your PATH environment variable that can be changed in system settings, or you may simply copy them into your System32 (x64) or SysWOW64 (x86) directories.

Run the Demo
--------------------
To run the demo of CreamyStrand, you may simply use the command line argument *-f [scene_file]* to specify the scene to be loaded. For example, you may type

./StrandSimulatorApp -f assets/unit_tests/drag/default_strand_cream.xml

to run the simulation of shaving cream poured onto seven pinned strands. 

All the parameters can be modified offline in the scene description XML files.

*NOTE*: For running on Windows, please set the *Working Directory* to be the directory containing the assets folder. Some examples may read data from the files in this folder.

USAGE:

   ./StrandSimulatorApp [-l <0-2>] [-d <integer>] [-c <integer>] [-o <integer>] -f <string> [--] [--version] [-h]

Where:

   -l <0-2>,  --statlog <0-2>
     Log runtime stats: 0: no stats, 1: timestep stats, 2: all stats

   -d <integer>,  --display <integer>
     Run the simulation with display enabled if 1, without if 0

   -c <integer>,  --checkpoint <integer>
     Between # steps several binary files are written to cache simulation state to

   -o <integer>,  --outputfile <integer>
     Between # steps several PLY files are written to save simulation state to

   -f <string>,  --file <string>
     (required)  XML file for a problem

   --,  --ignore_rest
     Ignores the rest of the labeled arguments following this flag.

   --version
     Displays version information and exits.

   -h,  --help
     Displays usage information and exits.

Surface Reconstruction and Rendering with Houdini
--------------------------------------------------------
Our simulator can generate PLY files that can be read back by SideFX Houdini.

You may run the demo with the "-o" option. For example, you may type

./StrandSimulatorApp -f assets/unit_tests/drag/default_strand_cream.xml -o 4

to generate data per 4 time steps (for this example we use 0.001s for the simulation time step, and 0.004s for the rendering time step. Hence 4 time steps is taken for the data generation). The simulation code will create a folder with the name of the scene file and timestamp, which contains all the generated data.

For more information and tutorials of Houdini, please visit the SideFX website (https://www.sidefx.com/).

Contact
-----------
Please contact the author (fyun@acm.org) for questions and bug report, or consulting for the usage of this code base.

BibTex Citation
----------------------
@article{fei2019mmc,  
 author = {Fei, Yun (Raymond) and Batty, Christopher and Grinspun, Eitan and Zheng, Changxi},  
 title = {A Multi-scale Model for Coupling Strands with Shear-Dependent Liquid},  
 journal = {ACM Trans. Graph.},  
 issue_date = {November 2019},  
 volume = {38},  
 number = {6},  
 month = nov,  
 year = {2019},  
 pages = {1:1--1:20},  
 articleno = {1},  
 numpages = {20},  
 url = {http://doi.acm.org/10.1145/},  
 doi = {10.1145/3355089.3356532},  
 acmid = {3356532},  
 publisher = {ACM},  
 address = {New York, NY, USA}  
}
