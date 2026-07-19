% Regression test for the spherical (r, theta, phi) 3D operator set and the
% coordinate-system class hierarchy. Run with:  run Tests/FVTool_spherical3D_test
% Written for the coordinate-system refactor. See the license file.
clc

%% Class hierarchy: correct subclass, integer dimension, coordsystem, tag
checks = { ...
  createMesh1D(5,1),                     'Mesh1D',            1, 'cartesian',   '1D'; ...
  createMeshCylindrical1D(5,1),          'MeshCylindrical1D', 1, 'cylindrical', 'Cylindrical1D'; ...
  createMeshSpherical1D(5,1),            'MeshSpherical1D',   1, 'spherical',   'Spherical1D'; ...
  createMesh2D(5,4,1,1),                 'Mesh2D',            2, 'cartesian',   '2D'; ...
  createMeshCylindrical2D(5,4,1,1),      'MeshCylindrical2D', 2, 'cylindrical', 'Cylindrical2D'; ...
  createMeshRadial2D(5,4,1,pi),          'MeshRadial2D',      2, 'radial',      'Radial2D'; ...
  createMesh3D(3,4,5,1,1,1),             'Mesh3D',            3, 'cartesian',   '3D'; ...
  createMeshCylindrical3D(3,4,5,1,pi,1), 'MeshCylindrical3D', 3, 'cylindrical', 'Cylindrical3D'; ...
  createMeshSpherical3D(3,4,5,1,pi,2*pi),'MeshSpherical3D',   3, 'spherical',   'Spherical3D'; ...
};
for i=1:size(checks,1)
  mi = checks{i,1};
  assert(strcmp(class(mi), checks{i,2}), 'wrong class for %s', checks{i,5});
  assert(isa(mi,'MeshStructure'), 'not a MeshStructure: %s', checks{i,5});
  assert(mi.dimension==checks{i,3}, 'wrong dimension for %s', checks{i,5});
  assert(strcmp(mi.coordsystem, checks{i,4}), 'wrong coordsystem for %s', checks{i,5});
  assert(strcmp(geometryTag(mi), checks{i,5}), 'wrong geometryTag for %s', checks{i,5});
end
disp('Mesh class hierarchy verified for all 9 geometries.');

%% Every Spherical3D operator must build without error
m = createMeshSpherical3D(6, 5, 7, 1.0, pi, 2*pi);
D = createCellVariable(m, 1.0); Df = harmonicMean(D);
u = createFaceVariable(m, [0.1 0.0 0.0]);
phi = createCellVariable(m, 1.0);
BC = createBC(m);
FL = @(r)(r+abs(r))./(1+abs(r)); % van Leer limiter
n = prod(m.dims+2);
assert(isequal(size(diffusionTerm(Df)), [n n]));
assert(isequal(size(convectionTerm(u)), [n n]));
assert(isequal(size(convectionUpwindTerm(u)), [n n]));
[Mtvd,~] = convectionTvdTerm(u, phi, FL); assert(isequal(size(Mtvd), [n n]));
divergenceTerm(u); gradientTerm(phi); gradientCellTerm(phi);
[Mbc,~] = boundaryCondition(BC); assert(isequal(size(Mbc), [n n]));
cellBoundary(phi.value(2:end-1,2:end-1,2:end-1), BC);
cellVolume(m);
disp('All Spherical3D operators built successfully.');

%% Physics: steady radial diffusion converges to analytic -1+2/r at 2nd order
errs = zeros(1,3); Nrs = [20 40 80];
for kk = 1:3
    Nr = Nrs(kk);
    mm = createMeshSpherical3D(linspace(1,2,Nr+1), linspace(0,pi,5), linspace(0,2*pi,5));
    Dfh = harmonicMean(createCellVariable(mm, 1.0));
    BCs = createBC(mm);
    BCs.left.a(:)=0;  BCs.left.b(:)=1;  BCs.left.c(:)=1;
    BCs.right.a(:)=0; BCs.right.b(:)=1; BCs.right.c(:)=0;
    [Mb, RHSb] = boundaryCondition(BCs);
    sol = solvePDE(mm, -diffusionTerm(Dfh) + Mb, RHSb);
    analytic = -1 + 2.0./mm.cellcenters.x;
    errs(kk) = max(abs(squeeze(sol.value(2:end-1,2,2)) - analytic))/max(abs(analytic));
end
fprintf('radial-diffusion errors: %.2e -> %.2e -> %.2e\n', errs);
assert(errs(3) < 2e-4, 'spherical3D diffusion not accurate enough');
assert(errs(1)/errs(2) > 3.5 && errs(2)/errs(3) > 3.5, 'not 2nd-order convergent');
disp('Spherical3D diffusion validated (2nd-order, matches analytic).');

disp('=== FVTool_spherical3D_test PASSED ===');
