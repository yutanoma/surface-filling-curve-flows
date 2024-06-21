# Overview

This document describes the data formats and directory structure of the example files for "Surface-Filling Curve Flows via Implicit Medial Axes" by Yuta Noma, Silvia Sellán, Nicholas Sharp, Karan Singh, and Alec Jacobson.

The directory structure was inspired by the official code release of ["Repulsive Curves"](https://github.com/chrisyu-cs/repulsive-curves) by Chris Yu, Henrik Schumacher, and Keenan Crane.

# Directory structure

Example files are contained in `models/surfaces`. Each example is contained
in a subdirectory with an identical structure:

```
models/surfaces/airplane:
| scene.txt
| plane_holes.obj
```

The file `scene.txt`, detailed below, contains the "scene" describing the example. The scene file and the `.obj` file MUST be in the same directory.

# Surfaces

Surface files are stored as Wavefrom OBJ files. In addition to the common triangle mesh, one can also add a line that describes the initial curve for the flow

```
l j1 j2 j3 j4 j5 j6
```

as well as a list of fixed points that will not move throughout the flow

```
fixed j1 j6
```

# Scene files

Each line of a scene files describes either the scene geometry or an initial parameter.

For more details, please refer to `modules/src/scene_file.cpp`.

## Mesh

Each scene file must specify one and only one surface via a line

`mesh mesh.obj`

where `mesh` is a fixed keyword and a path is given to the mesh OBJ file.
Make sure that the mesh is in the same directory as the scene file.

## Curve

If necessary, the user can specify the initial curve via a line

`curve curve.obj`

where `curve` is a fixed keyword and a path is given to the curve OBJ file.
If not specified, the edges on the face with id=1 will be added.

## Scalar

The user can specify a spatially-varying alpha value to generate curves of varying density (see Figure 23 in the paper) via a line

`scalar a.txt`

where `a.txt` containes the alpha value per vertex (see models/surfaces/bob/bob_grayscale.txt for details).

To enable this, one must also add a line

`varying_alpha`

## Radius

`radius v`

Specifies the spacing of the surface-filling curve (see the paper for details). By default, it is set to 30.

## h

`h val`

Specifies the h value in the paper that handles the resolution of the curve.

## r_max

`rmax val`

Specifies the r_max value in the paper that handles the speed of the curve evolution.

## p, q

```
p 2
q 2
```

Specifies the order of the energy (see Appendix A.3).

## Additional energies

`field_aligned val`

Specifies to add a field alignment energy with a weight `val`, and

`bilaplacian val`

does to add a biharmonic bending energy with a weight `val`. The vector field computation is determined based on the `computeSmoothestVertexDirectionField` function in geometry central.

## Geodesic medial axis

`geodesic_medial_axis`

Switches the medial axis computation to use the geodesic medial axis (see Sec. 3.2.2 in the paper).
