% Validation test: transient heat diffusion with a Gaussian (photothermal)
% source in a 1D cylindrical mesh vs. the analytic Sheldon(1982) solution.
% Ported from PyFVTool/tests/test_cylindrical1D_diffusion_source_photothermal.py.
%   run Tests/FVTool_cylindrical1D_photothermal_test
%
% A Gaussian laser beam heats an absorbing medium.  PDE (radial only):
%   cp*rho dT/dt = qdot(r) + k laplacian(T),
%   qdot(r) = alpha S(r) = alpha (2P/(pi w^2)) exp(-2 r^2 / w^2).
% Analytic solution (Sheldon, S.J. et al., Appl. Opt. 1982, 21, 1663):
%   dT(r,t) = (2 P alpha)/(pi cp rho w^2) *
%             INT_0^t [1/(1+2t'/tc)] exp(-2(r^2/w^2)/(1+2t'/tc)) dt'
% with tc = (w^2 cp rho)/(4k).  The integral is evaluated numerically.
clc

P = 1.0e-3;      % laser power [W]
w = 50e-6;       % beam waist [m]
k = 0.598;       % [W/m/K]
rho = 998.23;    % [kg/m^3]
cp = 4184.0;     % [J/kg/K]
alpha = 1.0;     % absorption coefficient (set to 1)
tc = (w^2*cp*rho)/(4*k);

Sbeam = @(r) (2*P)/(pi*w^2) * exp(-2*r.^2/w^2);

%% ---- finite-volume solution -------------------------------------------
Nr = 100;
Lr = 10.0*w;
T0 = 0.0;
deltat = 0.001;
m = createMeshCylindrical1D(Nr, Lr);

fv_T = createCellVariable(m, T0);                       % no-flux boundaries
fv_dotq = createCellVariable(m, alpha*Sbeam(m.cellcenters.x));
k_face = harmonicMean(createCellVariable(m, k));
transcoeff = createCellVariable(m, cp*rho);

Mdiff = diffusionTerm(k_face);
srcRHS = constantSourceTerm(fv_dotq);
[Mbc, RHSbc] = boundaryCondition(createBC(m));

t = 0.0;
for i = 1:100
    % assembly: (Mt - Mdiff)*T = RHSt + RHSbc + srcRHS  encodes
    %   cp*rho dT/dt = k laplacian(T) + qdot
    [Mt, RHSt] = transientTerm(fv_T, deltat, transcoeff);
    fv_T = solvePDE(m, Mt - Mdiff + Mbc, RHSt + RHSbc + srcRHS);
    t = t + deltat;
end

%% ---- benchmark vs analytic --------------------------------------------
r_num = m.cellcenters.x;
phi_num = fv_T.value(2:Nr+1);
% Sheldon(1982) analytic profile: integral over t' by fine trapezoidal rule
A = (2*P*alpha)/(pi*cp*rho*w^2);
tp = linspace(0, t, 4000)';                     % column: integration variable
r2w2 = (r_num(:).^2/w^2)';                       % row: one per radial position
denom = 1.0 + 2*tp/tc;                           % column
integrand = (1.0./denom) .* exp((-2*r2w2)./denom);   % [Nt x Nr]
phi_an = A*trapz(tp, integrand, 1);              % 1 x Nr
norm_err = (phi_num(:) - phi_an(:))/max(phi_an);
fprintf('cyl-1D photothermal: max norm. error = %.3e (t=%.3f, tc=%.3f ms)\n', ...
        max(abs(norm_err)), t, tc*1000);
assert(all(abs(norm_err) < 2.5e-3), 'cylindrical1D photothermal error too large');

disp('=== FVTool_cylindrical1D_photothermal_test PASSED ===');
