% Validation test: diffusion of an initial uniform sphere into an infinite
% medium (1D spherical) vs. Crank's analytic error-function solution.
% Ported from PyFVTool/tests/test_spherical1D_diffusion.py.
%   run Tests/FVTool_spherical1D_diffusion_test
%
% Reference: J. Crank, "The Mathematics of Diffusion", 2nd Ed. (1975),
%            pp. 29-30, eqn (3.8), Fig. 3.1.  Initial condition
%   c(r<=a)=C0, c(r>a)=0,  free diffusion, analytic:
%   C = C0/2 [erf((a-r)/2s) + erf((a+r)/2s)]
%       - C0/r * sqrt(Dt/pi) [exp(-(a-r)^2/4Dt) - exp(-(a+r)^2/4Dt)]
% with s = sqrt(Dt).
%
% NOTE: a fine mesh is required for good agreement (as in the original
% PyFVTool example); this test therefore takes a few seconds.
clc

a = 1.0;          % radius of the initial sphere
C0 = 1.0;         % initial concentration inside the sphere
D_val = 1.0;      % diffusion coefficient

C_sphere_infmed = @(r,t) 0.5*C0*(erf((a-r)./(2*sqrt(D_val*t))) + erf((a+r)./(2*sqrt(D_val*t)))) ...
    - (C0./r).*sqrt(D_val*t/pi).*(exp(-(a-r).^2/(4*D_val*t)) - exp(-(a+r).^2/(4*D_val*t)));

%% ---- finite-volume solution -------------------------------------------
R = 10.0;         % domain radius (models an infinite medium)
Nr = 2000;        % number of cells (fine grid needed for accuracy)
m = createMeshSpherical1D(Nr, R);

D = createCellVariable(m, D_val);
alfa = createCellVariable(m, 1.0);
Mdiff = diffusionTerm(harmonicMean(D));

r_fvm = m.cellcenters.x;
c = createCellVariable(m, 0.0);
c.value(2:Nr+1) = C0*(r_fvm < a);     % initial condition (no-flux boundaries)

BC = createBC(m);
[Mbc, RHSbc] = boundaryCondition(BC);

deltat = 0.0625/20;
t = 0.0;
% integrate to t = 1.0 (20 + 60 + 240 = 320 steps)
for nstep = [20 60 240]
    for n = 1:nstep
        [Mt, RHSt] = transientTerm(c, deltat, alfa);
        c = solvePDE(m, Mt - Mdiff + Mbc, RHSt + RHSbc);
        t = t + deltat;
    end
end

%% ---- benchmark vs analytic (at t = 1.0) --------------------------------
fvmdiff = c.value(2:Nr+1) - C_sphere_infmed(r_fvm, 1.0);
fprintf('sph-1D transient diffusion: max abs. diff = %.3e (t=%.3f)\n', ...
        max(abs(fvmdiff)), t);
assert(all(abs(fvmdiff) < 5e-4), 'spherical1D diffusion error too large');

disp('=== FVTool_spherical1D_diffusion_test PASSED ===');
