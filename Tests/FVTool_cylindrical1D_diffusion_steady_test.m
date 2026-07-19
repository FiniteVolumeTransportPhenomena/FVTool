% Validation test: steady heat diffusion with an internal source in a 1D
% cylindrical mesh vs. the analytic parabolic profile.
% Ported from PyFVTool/tests/test_cylindrical1D_diffusion_steady.py.
%   run Tests/FVTool_cylindrical1D_diffusion_steady_test
%
% PDE:   0 = k * laplacian(T) + S           (cylindrical, radial only)
% with Dirichlet outer wall T(R) = 0 and a symmetry (no-flux) centre,
% analytic solution:   T(r) = S/(4k) * (R^2 - r^2)
clc

R = 2.1;        % cylinder radius
k_val = 3.7;    % heat transfer coefficient
S_val = 4.2;    % volumetric source strength
T_outer = 0.0;  % outer-wall (Dirichlet) temperature
Nr = 50;

m = createMeshCylindrical1D(Nr, R);
k = createCellVariable(m, k_val);
S = createCellVariable(m, S_val);

BC = createBC(m);                                 % centre stays no-flux
BC.right.a(:) = 0.0; BC.right.b(:) = 1.0; BC.right.c(:) = T_outer;

k_face = harmonicMean(k);
Mdiff = diffusionTerm(k_face);
[Mbc, RHSbc] = boundaryCondition(BC);

% steady state: Mdiff*T = -constantSourceTerm(S)  (see poisson_example.m)
T = solvePDE(m, Mdiff + Mbc, RHSbc - constantSourceTerm(S));

r = m.cellcenters.x;
Tnum = T.value(2:Nr+1);
T_an = (S_val/(4*k_val))*(R^2 - r.^2);

norm_err = (Tnum(:) - T_an(:))/max(T_an);
fprintf('cyl-1D steady diffusion+source: max norm. error = %.3e\n', max(abs(norm_err)));
assert(all(abs(norm_err) < 1e-3), 'cylindrical1D steady diffusion error too large');

disp('=== FVTool_cylindrical1D_diffusion_steady_test PASSED ===');
