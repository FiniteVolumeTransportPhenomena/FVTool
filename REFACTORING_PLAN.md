# FVTool refactoring plan: coordinate systems as classes

Branch: `refactor`

## 1. Problem

`MeshStructure` stores a single `dimension` field that overloads **two** independent
concepts — the spatial dimension and the coordinate system — into one floating-point code:

| code  | dimension | coordinate system        | mesh creator              |
|-------|-----------|--------------------------|---------------------------|
| `1`   | 1         | Cartesian                | `createMesh1D`            |
| `1.5` | 1         | Cylindrical (radial)     | `createMeshCylindrical1D` |
| `1.8` | 1         | Spherical                | `createMeshSpherical1D`   |
| `2`   | 2         | Cartesian                | `createMesh2D`            |
| `2.5` | 2         | Cylindrical (axisym r,z) | `createMeshCylindrical2D` |
| `2.8` | 2         | Radial / polar (r,θ)     | `createMeshRadial2D`      |
| `3`   | 3         | Cartesian                | `createMesh3D`            |
| `3.2` | 3         | Cylindrical (r,θ,z)      | `createMeshCylindrical3D` |
| `3.5` | 3         | Spherical (r,θ,φ)        | `createMeshSpherical3D`   |

### Why this is bad practice
- **Float equality dispatch.** ~45 files branch with `switch d ... case 1.8` or
  `if d==2.8`. The codes `1.8`, `2.8`, `3.2`, `3.5` are **not exactly representable in
  binary floating point**; dispatch works today only because the identical literal is
  reused everywhere. Any arithmetic, rounding, or reformatting silently breaks routing.
- **Silent fall-through.** No dispatcher has an `otherwise`/`else`, so an unmatched code
  returns an *unset* output variable instead of raising an error.
- **Two concepts, one number.** Code that only needs the spatial dimension (e.g.
  `harmonicMean`) must still enumerate every geometry code, and the meaning of `2.8` is
  not discoverable from the value.

### Confirmed breakage
- **Spherical 3D (`3.5`) is entirely dead.** No dispatcher has a `case 3.5` and **no
  `*Spherical3D` operator exists**. `createMeshSpherical3D` produces empty matrices.
- **`createMeshTilted2D`** is tagged plain `2`; its tilted metric is never dispatched.
- **`Examples/Tutorial/tracer_tvd.m`** branches on `m.dimension==1 / ==2`, missing the
  `1.5/2.5/2.8` variants (works only by luck on Cartesian meshes).

## 2. Target design — subclass hierarchy

Mirror the PyFVTool structure (`Grid1D`, `CylindricalGrid1D`, …) while keeping FVTool's
functional MATLAB style and Octave compatibility.

```
MeshStructure                 (base; holds dims/cellsize/cellcenters/facecenters/corners/edges)
├── Mesh1D                     coordsystem='cartesian',   dimension=1
│   ├── MeshCylindrical1D      coordsystem='cylindrical', dimension=1
│   └── MeshSpherical1D        coordsystem='spherical',   dimension=1
├── Mesh2D                     coordsystem='cartesian',   dimension=2
│   ├── MeshCylindrical2D      coordsystem='cylindrical', dimension=2
│   └── MeshRadial2D           coordsystem='radial',      dimension=2
└── Mesh3D                     coordsystem='cartesian',   dimension=3
    ├── MeshCylindrical3D      coordsystem='cylindrical', dimension=3
    └── MeshSpherical3D        coordsystem='spherical',   dimension=3
```

Two clean fields replace the float code:
- `dimension` — **integer** `1|2|3` (true spatial dimension).
- `coordsystem` — string `'cartesian'|'cylindrical'|'radial'|'spherical'`.

### Dispatch mechanism
MATLAB single-dispatches on the *first* argument, which for term functions is a
`FaceVariable`/`CellVariable`, not the mesh — so full method polymorphism on the mesh
would require ~7×12 method files. Instead we keep the thin free-function dispatchers and
switch on a **canonical tag** computed by a method on the base class:

```matlab
function tag = geometryTag(m)   % method on MeshStructure
    coord = [upper(m.coordsystem(1)) m.coordsystem(2:end)];  % 'Cartesian', ...
    if strcmp(m.coordsystem,'cartesian')
        tag = sprintf('%dD', m.dimension);                    % '1D','2D','3D'
    else
        tag = sprintf('%s%dD', coord, m.dimension);           % 'Cylindrical2D', ...
    end
end
```

The tag deliberately reproduces the existing implementation-file suffixes
(`1D`, `Cylindrical2D`, `Radial2D`, `Spherical3D`, …). Dispatchers become:

```matlab
switch geometryTag(u.domain)
    case '1D',            M = diffusionTerm1D(u);
    case 'Cylindrical1D', M = diffusionTermCylindrical1D(u);
    ...
    otherwise
        error('FVTool:unsupportedGeometry', ...
              'diffusionTerm: no implementation for %s', geometryTag(u.domain));
end
```

- **Category A** functions (need only spatial dimension: `harmonicMean`, `*Mean`,
  `internalCells`, `reshapeCell`, `createFaceVariable`, `cellVolume`, `solvePDE`, …) simply
  read the integer `dimension` — `switch m.dimension` with `case 1/2/3`.
- **Category B** functions (geometry-specific: diffusion, convection family, divergence,
  gradient, boundaryCondition, cellBoundary, visualization) use `geometryTag`.

Some tags intentionally map to a shared Cartesian implementation (e.g. `Cylindrical2D`
boundary conditions reuse `boundaryCondition2D`); the switch preserves those groupings
explicitly instead of relying on magnitude ranges.

### Octave & backward compatibility
- classdef single-inheritance and `obj@Super(...)` constructor chaining work in Octave ≥4.
- **Breaking change:** `dimension` now returns an integer (`2`) instead of a float (`2.5`).
  This is intentional. In-repo examples/tests that compared against float codes are updated.
  External user scripts that relied on `dimension==2.5` must switch to
  `strcmp(m.coordsystem,'cylindrical')`. Documented in CLAUDE.md/README.

## 3. Work breakdown

### Phase 1 — Class hierarchy (`Classes/`, `MeshGeneration/`)
1. Add `coordsystem` property + `geometryTag` method to `@MeshStructure`; keep base
   constructor accepting the field list, default `coordsystem='cartesian'`.
2. Add subclass files `Classes/Mesh1D.m`, `Mesh2D.m`, `Mesh3D.m`,
   `MeshCylindrical1D.m`, `MeshSpherical1D.m`, `MeshCylindrical2D.m`, `MeshRadial2D.m`,
   `MeshCylindrical3D.m`, `MeshSpherical3D.m` (each sets `dimension` + `coordsystem`).
3. Rewrite the 13 `createMesh*.m` to construct the right subclass with an integer
   dimension (drop every `MS.dimension=<float>` line). `createMeshTilted2D` returns a
   `Mesh2D` (documented as Cartesian until a tilted operator set exists).

### Phase 2 — Dispatcher refactor (~45 files)
Mechanical, file-by-file, validated by the existing test scripts:
- Category B → `switch geometryTag(...)` + `otherwise error`.
- Category A → `switch m.dimension` (`1|2|3`).
Covers: `Discretization/` (diffusion, convection×4, source, transient), `Calculus/`
(divergence, gradient×3, ddt), `Boundary/` (boundaryCondition, cellBoundary, createBC,
combineBC), `Utilities/` (all `*Mean`, internalCells, reshapeCell/InternalCell,
excludeGhostRHS, BC2GhostCells, maskCells, normCellVector, cellLocations, faceLocations),
`Solvers/` (solvePDE, solveExplicitPDE), `Visualization/` (visualizeCells,
visualizeCellVectors, visualizeMesh).

### Phase 3 — Spherical 3D operators (new files)
Port from PyFVTool (`diffusion.py`, `advection.py`, `boundary.py` `*Spherical3D`) using
the existing `Cylindrical3D` files as the MATLAB template:
- `Discretization/diffusionTermSpherical3D.m`
- `Discretization/convectionTermSpherical3D.m`, `convectionUpwindTermSpherical3D.m`,
  `convectionTvdTermSpherical3D.m`, `convectionTvdRHSSpherical3D.m`
- `Calculus/divergenceTermSpherical3D.m`, `gradientTermSpherical3D.m`
- `Boundary/boundaryConditionSpherical3D.m`, `cellBoundarySpherical3D.m`
- `Visualization/visualizeCellsSpherical3D.m` (+ register in visualize dispatchers)
- Add the `Spherical3D` case to every Phase-2 dispatcher.

### Phase 4 — Examples & tests
- Fix `tracer_tvd.m` and any other example branching on float codes.
- Add `Examples/Tutorial/diffusion_spherical3D.m`.
- Extend `Tests/FVTool_functions_uniform_test.m` / `_nonuniform_test.m` to construct and
  exercise every mesh subclass (including Spherical3D).

### Phase 5 — Validation & docs
- Run the three test scripts in MATLAB (`/usr/local/bin/matlab -batch`).
- Validate Spherical3D diffusion against an analytic solution (steady radial diffusion in a
  sphere reduces to the known 1/r profile; also cross-check against `Spherical1D`).
- Update `CLAUDE.md` (dimension-dispatch section) and `README.md`.

## 4. Risks
- **Octave classdef inheritance** edge cases — mitigated by keeping subclasses thin
  (properties + trivial constructor only; all logic stays in free functions).
- **Spherical3D numerics** — highest risk; mitigated by porting a validated reference and
  checking against analytic + 1D-spherical results before declaring done.
- **Silent behavior changes** from the integer-`dimension` break — mitigated by the loud
  `otherwise` errors and the full test sweep.
