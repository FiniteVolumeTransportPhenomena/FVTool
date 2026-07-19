# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

FVTool is a finite-volume PDE toolbox for MATLAB/Octave. It discretizes and solves the general transient
convection-diffusion equation

```
alpha * d(phi)/dt + div(u*phi) + div(-D grad phi) + beta*phi + gamma = 0
```

with Dirichlet/Neumann/Robin/periodic boundary conditions, on 1D/2D/3D Cartesian, cylindrical, radial, and
spherical grids. Pure MATLAB `.m` files — no build step, no package manager, no compiled artifacts.

## Running things

There is no build or lint. Everything runs inside MATLAB or Octave (>= 4.0, needs `classdef`).

- **Set up the path (do this first, every session):** run `FVToolStartUp` from the repo root. It `addpath`s all
  code directories and optionally detects the AGMG/Factorize solvers and the PVT toolbox if present.
- **Run the test suite:** `Tests/FVTool_functions_uniform_test.m` and `Tests/FVTool_functions_nonuniform_test.m`
  are scripts that exercise the discretization functions across all mesh types; `Tests/readme_test.m` runs the
  README example. Run a test by executing the script (e.g. `run Tests/FVTool_functions_uniform_test`). There is
  no per-function test runner — tests are scripts, not a framework.
- **Learn the API by example:** `Examples/Tutorial/` has focused single-equation demos; `Examples/Advanced/`
  has coupled/nonlinear reservoir-engineering solvers. `FVTdemo.m` is a large annotated walkthrough.

## Architecture

### The solve pattern
Every solution follows the same assembly pipeline. Understand this and the whole library follows:

1. **Create a mesh** — `createMesh2D(Nx,Ny,Lx,Ly)`, `createMeshCylindrical3D(...)`, etc. Returns a
   `MeshStructure`.
2. **Create variables on the mesh** — `createCellVariable` (values at cell centers, includes ghost cells),
   `createFaceVariable` / `createCellVector` (values at faces).
3. **Set boundary conditions** — `createBC(mesh)` returns a `BoundaryCondition` with per-boundary `a`, `b`, `c`
   arrays (default all-Neumann/zero-flux); edit `BC.left.a`, `BC.right.b`, etc.
4. **Assemble term matrices** — each PDE term is a function returning a sparse coefficient matrix `M` (and an
   RHS where relevant): `diffusionTerm`, `convectionTerm`/`convectionUpwindTerm`/`convectionTvdTerm`,
   `transientTerm`, `linearSourceTerm`, `constantSourceTerm`, and `boundaryCondition(BC)` -> `[Mbc, RHSbc]`.
5. **Sum matrices and RHS**, then `phi = solvePDE(mesh, M, RHS)`. `solvePDE` does `M\RHS` by default and
   reshapes the flat solution back into a `CellVariable` (with ghost cells). `solveExplicitPDE` for explicit steps.
6. **Visualize** — `visualizeCells(phi)`, `visualizeCellVectors(...)`.

Coefficient functions take face values, so cell-centered coefficients must first be interpolated to faces with
`harmonicMean` / `arithmeticMean` / `geometricMean` / `linearMean` / `upwindMean` / `tvdMean`.

### The dimension-dispatch convention (most important internal pattern)
User-facing functions are thin dispatchers. `diffusionTerm.m`, `convectionTerm.m`, `divergenceTerm.m`,
`boundaryCondition.m`, etc. read `D.domain.dimension` and `switch` to a geometry-specific implementation.
The dimension is encoded as a **fractional code**, not just the integer dimension count:

- `1` = 1D Cartesian, `1.5` = cylindrical/axisymmetric 1D, `1.8` = spherical 1D
- `2` = 2D Cartesian, `2.5` = cylindrical 2D, `2.8` = radial 2D (r, theta)
- `3` = 3D Cartesian, `3.2` = cylindrical 3D

So `diffusionTermCylindrical2D.m`, `diffusionTermRadial2D.m`, `diffusionTermSpherical1D.m` are the real
implementations; `diffusionTerm.m` just routes to them. **When adding or fixing a discretization, edit the
geometry-specific `*1D/2D/3D/Cylindrical*/Radial*/Spherical*` file and make sure the dispatcher's `switch`
covers the code.** The set of implemented geometries varies per operator — check which specific files exist
before assuming a geometry is supported.

### Operator overloading
`CellVariable`, `FaceVariable`, and `CellVector` are `classdef` types under `Classes/@.../`. Arithmetic and
comparison operators (`plus`, `minus`, `times`, `rdivide`, `power`, `gt`, `le`, ...) are overloaded per class
so 1D-shaped user code carries over unchanged to 2D/3D. The overloads operate on the `.value` field and return
the same class, preserving `.domain`. This is why example scripts can write `D.*grad(phi)` naturally regardless
of dimension.

### Directory map
- `MeshGeneration/` — `createMesh*`, `createCellVariable`, `createFaceVariable`, `cellVolume`
- `Classes/` — the four OO types and their operator overloads
- `Discretization/` — the PDE term coefficient matrices (the dispatch + geometry files above)
- `Boundary/` — `createBC`, `boundaryCondition`, ghost-cell handling (`cellBoundary`)
- `Calculus/` — `divergenceTerm`, `gradientTerm`, `ddtTerm` (post-processing differential operators)
- `Solvers/` — `solvePDE`, `solveExplicitPDE` (drop `AGMG_*` or `Factorize` dirs here to enable them)
- `Utilities/` — face-interpolation means, `fluxLimiter`, `RhieChow`, reshaping helpers
- `Visualization/` — `visualizeCells*`, `visualizeCellVectors*`, `visualizeMesh*`
- `PhysicalProperties/`, `FieldGeology/`, `Physics/` — domain-specific correlations (CO2/brine PVT, relative
  permeability, capillary pressure, Buckley-Leverett analytics). Optional; `FVToolStartUp` tolerates them being
  absent.

## Conventions when editing

- Follow the existing file-per-function, dispatcher + geometry-variant layout. Adding a scheme means touching
  the dispatcher and adding parallel files for each supported geometry.
- Ghost cells: `CellVariable.value` includes a ghost layer (dims are `N+2`). Interior/boundary handling lives
  in `Boundary/` and `Utilities/internalCells.m`, `excludeGhostRHS.m`, `BC2GhostCells.m` — reuse these rather
  than reindexing by hand.
- Keep changes Octave-compatible (no MATLAB-only toolbox dependencies in core code).
- Function headers use the SYNOPSIS/PARAMETERS/RETURNS/EXAMPLE/SEE ALSO doc block; match it for new functions.
