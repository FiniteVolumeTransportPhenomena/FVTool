% Steady-state diffusion in a spherical shell (r, theta, phi), solved on a
% 3D spherical mesh. With Dirichlet values on the inner and outer radii and
% zero-flux (default Neumann) elsewhere, the solution is purely radial and
% matches the analytic profile
%
%     phi(r) = c1 + c2/r,
%
% so that for phi(r1)=phiL and phi(r2)=phiR:
%     c2 = (phiL-phiR)/(1/r1-1/r2),   c1 = phiL - c2/r1.
%
% This example exercises the MeshSpherical3D operator set.
%
% Written for the coordinate-system refactor. See the license file.

clc; clear;

% geometry of the spherical shell
r1 = 1.0;      % inner radius
r2 = 2.0;      % outer radius
phiL = 1.0;    % value at inner radius
phiR = 0.0;    % value at outer radius

Nr = 40; Ntheta = 8; Nphi = 8;
rfaces   = linspace(r1, r2, Nr+1);
thfaces  = linspace(0, pi, Ntheta+1);
phifaces = linspace(0, 2*pi, Nphi+1);
m = createMeshSpherical3D(rfaces, thfaces, phifaces);

% unit diffusion coefficient
D = createCellVariable(m, 1.0);
Dface = harmonicMean(D);

% boundary conditions: Dirichlet on r (left = inner, right = outer),
% zero-flux (default) on theta and phi
BC = createBC(m);
BC.left.a(:)  = 0; BC.left.b(:)  = 1; BC.left.c(:)  = phiL;
BC.right.a(:) = 0; BC.right.b(:) = 1; BC.right.c(:) = phiR;

% assemble and solve  div(-D grad phi) = 0
Mdiff = diffusionTerm(Dface);
[Mbc, RHSbc] = boundaryCondition(BC);
c = solvePDE(m, -Mdiff + Mbc, RHSbc);

% compare with the analytic radial profile
rp = m.cellcenters.x;
c2 = (phiL - phiR)/(1/r1 - 1/r2);
c1 = phiL - c2/r1;
analytic = c1 + c2./rp;
numeric  = squeeze(c.value(2:end-1, 2, 2));   % any (theta,phi) slice
fprintf('max relative error vs analytic phi(r)=c1+c2/r : %.3e\n', ...
    max(abs(numeric - analytic))/max(abs(analytic)));

% visualize
figure; visualizeCells(c); title('Spherical 3D diffusion');
