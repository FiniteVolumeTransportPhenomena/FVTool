% Validation test: the Mason-Weaver equation (sedimentation + diffusion in a
% closed 1D column) vs. analytic expectations for mass conservation and the
% steady-state amplitude.
% Ported from PyFVTool/tests/test_pdesolver_mason_weaver.py.
%   run Tests/FVTool_masonWeaver_test
%
% Reference: Midelet et al., Part. Part. Syst. Charact. 2017, 34, 1700095.
% PDE:  dc/dt = D d2c/dz2 + sg dc/dz,  closed (no-flux) ends.
% Checks:
%   (1) total amount of c is conserved (closed system);
%   (2) the solution approaches the analytic steady-state amplitude
%       B = z_max / ( z0 (1 - exp(-z_max/z0)) ),  z0 = D/sg.
clc

z_max = 1.0;
D_coeff = 0.015;
sg = 0.2;
Nx = 100;
dt = 0.01;
t_sim = 10.0;
maxdev_ppq = 1000.0;   % max allowed rel. deviation in parts per 1e15

m = createMesh1D(Nx, z_max);
BC = createBC(m);                       % all no-flux (closed column)
[Mbc, RHSbc] = boundaryCondition(BC);

c = createCellVariable(m, 1.0);

% advection field (sedimentation), with closed ends (no flow at extremities)
u = createFaceVariable(m, sg);
u.xvalue(1)   = 0.0;
u.xvalue(end) = 0.0;

% diffusion field
D = createFaceVariable(m, D_coeff);

Mdiff = diffusionTerm(D);
Mconv = convectionTerm(u);

total_c = domainInt(c);
total0 = total_c;
max_dev = 0.0;

it = 0;
while it*dt < t_sim
    [Mt, RHSt] = transientTerm(c, dt, 1.0);
    c = solvePDE(m, Mt - Mdiff + Mconv + Mbc, RHSt + RHSbc);
    it = it + 1;
    dev = 1e15*(domainInt(c) - total0)/total0;
    max_dev = max(max_dev, abs(dev));
end

% steady-state amplitude
z0 = D_coeff/sg;
B = z_max/(z0*(1.0 - exp(-z_max/z0)));

fprintf('Mason-Weaver: max mass deviation = %.3e ppq, max(c) = %.4f, 0.9*B = %.4f\n', ...
        max_dev, max(c.value), 0.9*B);
assert(max_dev < maxdev_ppq, 'Mason-Weaver mass not conserved');
assert(max(c.value) > 0.9*B, 'Mason-Weaver did not reach steady-state amplitude');

disp('=== FVTool_masonWeaver_test PASSED ===');
