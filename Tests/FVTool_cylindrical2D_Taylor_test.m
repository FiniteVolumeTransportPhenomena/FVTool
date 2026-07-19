% Validation test: convective (Taylor) dispersion of a solute slug in
% Poiseuille flow through a cylindrical tube (2D cylindrical), compared with
% G.I. Taylor's analytic radially-averaged profile.
% Ported from PyFVTool/tests/test_cylindrical2D_convection_Taylor.py.
%   run Tests/FVTool_cylindrical2D_Taylor_test
%
% Reference: G.I. Taylor, Proc. Roy. Soc. A 1953, 219, 186 (eqn "A3").
% Purely convective transport of an initial slug; the radially integrated
% concentration profile along the tube is compared to the analytic result.
%
% NOTE: 40x500 grid, 1000 time steps (a few seconds to run).
clc

Lr = 7.5e-05;            % [m] tube radius
Lz = 0.3;               % [m] tube length
umax = 2*9.4314e-3;     % [m/s] peak (= 2x mean) Poiseuille velocity
Nr = 40;  Nz = 500;
loadix0 = 20; loadix1 = 40;   % initial-slug column range (0-based interior)
deltat = 0.01;

m = createMeshCylindrical2D(Nr, Nz, Lr, Lz);

% Poiseuille velocity field: only axial (z) component, parabolic in r
r = m.cellcenters.x;
uu = umax*(1 - (r.^2)/Lr^2);
u = createFaceVariable(m, [0 0]);
u.yvalue = repmat(uu(:), 1, Nz+1);      % axial face velocities (Nr x Nz+1)

% solution variable, default no-flux BCs
phi = createCellVariable(m, 0.0);
% initial condition: unit slug over interior columns loadix0..loadix1-1
% (0-based interior -> +2 for FVTool full ghosted array)
phi.value(2:Nr+1, (loadix0+2):(loadix1+1)) = 1.0;

Mconv = convectionUpwindTerm(u);
[Mbc, RHSbc] = boundaryCondition(createBC(m));
v_int = internalCells(cellVolume(m));    % Nr x Nz interior cell volumes

% radial integral of the concentration for every axial position -> 1 x Nz
integral_dr = @(p) sum(v_int .* p.value(2:Nr+1, 2:Nz+1), 1);

% record the initial radial-integral profile (for the analytic amplitude C_0)
phiprof0 = integral_dr(phi);

% ---- time stepping (200 + 300 + 500 = 1000 steps) ----------------------
t = 0.0;
for nstep = [200 300 500]
    for i = 1:nstep
        [Mt, RHSt] = transientTerm(phi, deltat, 1.0);
        phi = solvePDE(m, Mt + Mconv + Mbc, RHSt + RHSbc);
        t = t + deltat;
    end
end
phiprof = integral_dr(phi);

%% ---- analytic Taylor A3 profile ---------------------------------------
zf = m.facecenters.y;                    % axial face centers (Nz+1)
DX = zf(loadix0+1);                       % 0-based -> 1-based
X  = zf(loadix1+1) - DX;
C_0 = phiprof0((loadix0+loadix1)/2 + 1);  % 0-based interior idx 30 -> 1-based 31
z_num = m.cellcenters.y;                  % axial cell centers (Nz)

xs = z_num - DX;
c_an = zeros(size(xs));
assert(t >= X/umax, 'analytic branch t < X/u0 not implemented');
for ix = 1:numel(xs)
    x = xs(ix);
    if (x >= 0) && (x < X)
        c_an(ix) = C_0 * x/(umax*t);
    elseif (x >= X) && (x < umax*t)
        c_an(ix) = C_0 * X/(umax*t);
    elseif (x >= umax*t) && (x < umax*t + X)
        c_an(ix) = C_0 * (X + umax*t - x)/(umax*t);
    else
        c_an(ix) = 0.0;
    end
end

norm_err = (c_an(:) - phiprof(:))/max(c_an);
% benchmark over the mid-range of the tube (Python slice [Nz/3 : Nz/2])
idx = (floor(Nz/3)+1):floor(Nz/2);
fprintf('cyl-2D Taylor dispersion: max norm. error (mid-tube) = %.3e (t=%.2f s)\n', ...
        max(abs(norm_err(idx))), t);
assert(all(abs(norm_err(idx)) < 1.5e-3), 'cylindrical2D Taylor dispersion error too large');

disp('=== FVTool_cylindrical2D_Taylor_test PASSED ===');
