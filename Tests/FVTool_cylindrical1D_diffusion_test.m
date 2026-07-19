% Validation test: transient heat/mass diffusion into a solid cylinder from
% a step change at the outer surface, vs. Crank's analytic Bessel series.
% Ported from PyFVTool/tests/test_cylindrical1D_diffusion.py.
%   run Tests/FVTool_cylindrical1D_diffusion_test
%
% Reference: J. Crank, "The Mathematics of Diffusion", 2nd Ed. (1975),
%            p. 78, eqn (5.22).  A cylinder of radius a, initially at C=0,
%            surface held at C=1 for t>0:
%   C(r,t) = 1 - (2/a) * sum_n exp(-D alpha_n^2 t) J0(r alpha_n)
%                                   / ( alpha_n J1(a alpha_n) )
% where a*alpha_n are the positive zeros of J0.
clc

a_val = 2.9;      % cylinder radius
D_val = 1.9;      % diffusion coefficient
Nterm = 30;       % number of series terms
c_outer = 1.0;    % surface concentration
c0 = 0.0;         % initial concentration

% ---- zeros of J0 (McMahon asymptotic guess refined by Newton's method) --
% J0'(x) = -J1(x), so Newton step is x <- x + J0(x)/(-J1(x)) ... solving
% J0(x)=0:  x <- x - J0(x)/(-J1(x)) = x + J0(x)/J1(x).
aalp = zeros(Nterm,1);
for n = 1:Nterm
    x = (n - 0.25)*pi;                 % McMahon initial estimate
    for it = 1:60
        x = x + besselj(0,x)/besselj(1,x);
    end
    aalp(n) = x;
end
alpha = aalp/a_val;

crank522 = @(r,t) 1.0 - (2.0/a_val) * ...
    sum( exp(-D_val*alpha.^2*t) .* besselj(0, r*alpha) ...
         ./ (alpha .* besselj(1, aalp)) );

%% ---- finite-volume solution -------------------------------------------
Nr = 50;
deltat = 0.001;
m = createMeshCylindrical1D(Nr, a_val);
c = createCellVariable(m, c0);

BC = createBC(m);                                  % centre: no-flux
BC.right.a(:) = 0.0; BC.right.b(:) = 1.0; BC.right.c(:) = c_outer;
[Mbc, RHSbc] = boundaryCondition(BC);

D = createCellVariable(m, D_val);
alfa = createCellVariable(m, 1.0);
Mdiff = diffusionTerm(harmonicMean(D));

t = 0.0;
for i = 1:1000
    [Mt, RHSt] = transientTerm(c, deltat, alfa);
    c = solvePDE(m, Mt - Mdiff + Mbc, RHSt + RHSbc);
    t = t + deltat;
end

%% ---- benchmark vs analytic --------------------------------------------
r_num = m.cellcenters.x;
phi_num = c.value(2:Nr+1);
phi_an = arrayfun(@(rr) crank522(rr, t), r_num);
norm_err = (phi_num(:) - phi_an(:))/max(phi_an);
fprintf('cyl-1D transient diffusion: max norm. error = %.3e (t=%.3f)\n', ...
        max(abs(norm_err)), t);
assert(all(abs(norm_err) < 5e-4), 'cylindrical1D transient diffusion error too large');

disp('=== FVTool_cylindrical1D_diffusion_test PASSED ===');
