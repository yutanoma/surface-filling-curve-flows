## Surface-Filling Curve Flows via Implicit Medial Axes

Public code release for [Surface-Filling Curve Flows via Implicit Medial Axes](https://www.dgp.toronto.edu/projects/surface-filling-curves/).

Surface-Filling Curve Flows via Implicit Medial Axes

Yuta Noma, Silvia Sellán, Nicholas Sharp, Karan Singh, Alec Jacobson

ACM Transaction on Graphics (Proceedings of SIGGRAPH North America 2024)

[Project page](https://www.dgp.toronto.edu/projects/surface-filling-curves/)

## Requirements

- Eigen: 3.4.0

## Build the code

First, install all the dependencies with `git submodule update --init --recursive`.

**Unix-like machines**: configure (with cmake) and compile
```
cd /path/to/directory
mkdir build
cd build
cmake ..
make -j6
```

**Windows / Visual Studio**

Install CMake, and use either the CMake GUI or the command line interface (as on unix) to generate a Visual Studio solution.  Build the solution with Visual Studio.

## Run the code
```
./bin/curve_on_surface ../models/surfaces/scorpion/scene.txt
```
By pressing the space bar or checking "run loop", it will produce a surface-filling curve!
