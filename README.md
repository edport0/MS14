# MS14 Project - 2D Mesh Generation and Adaptation

This repository contains a small C codebase used for the MS14 project on 2D triangular meshes.
The current work focuses on:
- naive point insertion with optional quality-based edge swaps,
- Delaunay point insertion with a cavity operator,
- mesh quality analysis and simple benchmarking helpers.

## Layout

- `src/`: C sources, build targets, and helper scripts
- `data/`: input meshes used for tests and comparisons
- `bin/`: external binaries provided with the course

## Build

Builds are driven by [`src/Makefile`](src/Makefile).

Compiler settings:
- compiler: `gcc`
- compile options: `-O2 -Wall`
- include path: `-I.`
- link option: `-lm`

Run the commands from `src/`.

Main targets:

- `make main_mesh`
  Builds `main_mesh`
- `make main_quality`
  Builds `main_quality`
- `make cavity`
  Builds `cavity`
- `make triangulation_naive`
  Builds `triangulation_naive`
- `make triangulation_delaunay`
  Builds `triangulation_delaunay`
- `make swap_test`
  Builds `swap_test`
- `make clean`
  Removes generated executables and object files

## Main Programs

### `main_mesh.c` -> `main_mesh`

General mesh inspection tool.

What it does:
- reads a mesh,
- builds neighbors with the quadratic and hash-table methods,
- computes triangle quality,
- writes:
  - `quality.solb`
  - `metric.solb`


```powershell
./main_mesh ..\data\square_unity_test.mesh
```

### `main_quality.c` -> `main_quality`

Quality-only post-processing tool.

What it does:
- reads a mesh,
- computes triangle quality,
- prints Q2 min / mean / max,
- writes:
  - `quality.solb`



```powershell
./main_quality triangulation_naive_out.mesh
./main_quality triangulation_delaunay_out.mesh
```

### `main_triangulation_naive.c` -> `triangulation_naive`

Naive incremental point insertion in the unit square.

What it does:
- inserts random points,
- localizes points either with:
  - `walk` mode using neighbor walk,
  - `full` mode using full mesh scan,
- splits triangles or edges,
- applies a local quality-based edge-swap improvement,
- writes `triangulation_naive_out.mesh`.

Usage:

```powershell
./triangulation_naive <mesh> [npts] [walk|full]
```

Examples:

```powershell
./triangulation_naive ..\data\square_unity_test.mesh 200 walk
./triangulation_naive ..\data\square_unity_test.mesh 200 full
```

Timing output:
- prints a human-readable insertion time,
- prints one compact line:
  - `TIMING naive walk <npts> <time>`
  - `TIMING naive full <npts> <time>`

### `main_triangulation_delaunay.c` -> `triangulation_delaunay`

Delaunay insertion based on the cavity operator.

What it does:
- inserts random points,
- localizes points,
- builds the Delaunay cavity,
- deletes cavity triangles,
- stars the cavity around the new point,
- updates neighbors locally,
- writes `triangulation_delaunay_out.mesh`.

Usage:

```powershell
./triangulation_delaunay <mesh> [npts]
```

Example:

```powershell
./triangulation_delaunay ..\data\square_unity_test.mesh 200
```

Timing output:
- prints a human-readable insertion time,
- prints one compact line:
  - `TIMING delaunay <npts> <time>`

### `main_cavity_test.c` -> `cavity`

Small debugging program for the cavity operator.

What it does:
- inserts one test point,
- prints cavity triangles,
- prints cavity boundary edges (`PilEdg`),
- applies cavity deletion and starring,
- prints the local retriangulation.

Usage:

```powershell
./cavity <mesh> [x y]
```

### `main_swap_test.c` -> `swap_test`

Small debugging program for edge swaps.

What it does:
- reads a mesh,
- builds neighbors,
- swaps one chosen interior edge,
- writes `swap_test_out.mesh`.

Usage:

```powershell
./swap_test <mesh> [iTri] [iEdg]
```

## Core Sources

### `mesh.c` / `mesh.h`

Main mesh data structures and algorithms:
- mesh reading/writing,
- bounding box and neighbors,
- hash table for adjacency,
- triangle quality functions,
- point location,
- naive split operators,
- cavity operator,
- local Delaunay insertion,
- edge swap and quality-based swap criterion.

### `libmesh6.c` / `libmesh6.h`

Mesh I/O support library for `.mesh` / `.sol` style files.

### `eigen.c`

Auxiliary numerical routines used by the project.



## Typical Workflow

1. Build the executable you need with `make`.
2. Generate a mesh:
   - `triangulation_naive`
   - or `triangulation_delaunay`
3. Inspect mesh quality with:
   - `main_mesh` for the full analysis path
   - `main_quality` for quality-only export
4. Use ViZiR to display:
   - the `.mesh`
   - the generated `.solb` quality field
