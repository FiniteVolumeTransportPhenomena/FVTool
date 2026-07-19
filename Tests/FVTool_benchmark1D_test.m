% Validation test: 1D Cartesian finite-volume solutions vs. analytic solutions.
% Ported from PyFVTool/tests/test_benchmark_1d.py.  Run with:
%   run Tests/FVTool_benchmark1D_test
%
% Three benchmarks:
%   1. Transient conduction, Dirichlet surface  -> error-function profile
%   2. Transient conduction, Neumann (fixed-flux) surface -> analytic profile
%   3. Steady convection-diffusion              -> exponential profile
%
% Each compares the numerical result against the closed-form solution and
% asserts a tight tolerance on the mean relative error.
clc

%% ---- 1 & 2: transient conduction in a semi-infinite slab ----------------
% PDE:  rho*c dT/dt = k d2T/dx2   (constant properties)
L  = 1.0;        % [m] domain length
k  = 20.0;       % [W/m/K]
rho = 8000.0;    % [kg/m^3]
cp = 500.0;      % [J/kg/K]
alfa = k/(rho*cp);   % thermal diffusivity
T0 = 300.0;      % [K] initial (and far-field) temperature
Ts = 350.0;      % [K] Dirichlet surface temperature
qs = 1000.0;     % [W/m^2] Neumann surface heat flux
t_sim = L^2/(20*alfa);
Nt = 50;
dt = t_sim/Nt;
Nx = 50;

% analytic profiles (semi-infinite medium)
T_an_dirichlet = @(x,t) (T0-Ts)*erf(x./sqrt(4*alfa*t))+Ts;
T_an_neumann   = @(x,t) T0 + qs/k*sqrt(4*alfa*t/pi).*exp(-x.^2/(4*alfa*t)) ...
                        - qs/k*x.*(1-erf(x./sqrt(4*alfa*t)));

for bc = {'Dirichlet','Neumann'}
    left_bc = bc{1};
    m = createMesh1D(Nx, L);
    BC = createBC(m);
    if strcmp(left_bc,'Dirichlet')
        BC.left.a(:) = 0.0;  BC.left.b(:) = 1.0;  BC.left.c(:) = Ts;
        T_analytic = @(x) T_an_dirichlet(x, t_sim);
    else
        % FVTool left BC reads  a*dphi/dx + b*phi = c  (gradient in +x).
        % Analytic dT/dx|0 = -qs/k, so with a=k we need c = -qs.
        BC.left.a(:) = k;    BC.left.b(:) = 0.0;  BC.left.c(:) = -qs;
        T_analytic = @(x) T_an_neumann(x, t_sim);
    end

    T = createCellVariable(m, T0, BC);            % initial condition
    alfa_face = harmonicMean(createCellVariable(m, alfa));
    Mdiff = diffusionTerm(alfa_face);
    [Mbc, RHSbc] = boundaryCondition(BC);

    t = 0.0;
    while t < t_sim - 1e-12
        t = t + dt;
        [Mt, RHSt] = transientTerm(T, dt, 1.0);
        T = solvePDE(m, Mt - Mdiff + Mbc, RHSt + RHSbc);
    end

    x = m.facecenters.x;
    Tface = linearMean(T);
    Tnum = Tface.xvalue;
    Tan  = T_analytic(x);
    er = sum(abs(Tnum(:)-Tan(:))./Tan(:))/Nx;
    fprintf('1D %-10s conduction: mean rel. error = %.3e\n', left_bc, er);
    assert(er <= 1e-3, '1D %s conduction error too large', left_bc);
end

%% ---- 3: steady convection-diffusion ------------------------------------
% u dc/dx = D d2c/dx2, analytic c = (1-exp(u x/D))/(1-exp(u L/D))
L  = 1.0; Nx = 50;
m = createMesh1D(Nx, L);
BC = createBC(m);
BC.left.a(:)  = 0; BC.left.b(:)  = 1; BC.left.c(:)  = 0;
BC.right.a(:) = 0; BC.right.b(:) = 1; BC.right.c(:) = 1;
x = m.cellcenters.x;

D_val = -1;  u_val = -10.0;
D_face = harmonicMean(createCellVariable(m, D_val));
u_face = createFaceVariable(m, u_val);

Mconv = convectionTerm(u_face);
Mdiff = diffusionTerm(D_face);
[Mbc, RHSbc] = boundaryCondition(BC);
c = solvePDE(m, Mconv - Mdiff + Mbc, RHSbc);

c_an = (1-exp(u_val*x/D_val))./(1-exp(u_val*L/D_val));
er = sum(abs(c_an(:)-c.value(2:Nx+1)))/Nx;
fprintf('1D steady convection-diffusion: mean abs. error = %.3e\n', er);
assert(er <= 1e-3, '1D convection error too large');

disp('=== FVTool_benchmark1D_test PASSED ===');
