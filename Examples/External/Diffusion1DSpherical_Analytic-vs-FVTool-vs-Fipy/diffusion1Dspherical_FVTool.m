%% Transient diffusion equation evaluated using FVTool
%
% adapted from A.A. Eftekhari by M.H.V. Werts (2020) 
%
% GNU Octave 4.2.2 was used for development
%
% This script is intended to be run from the command line interface. 
% It is called in this manner from within the accompanying Jupyter 
% Python notebook. 
%
% It should be inside the 'FVTool' directory tree as downloaded/cloned
% from Github, under 
%   './Examples/External/Diffusion1DSpherical_Analytic-vs-FVTool-vs-Fipy'
%
% Script calculates diffusion in a 1D spherical geometry for an 'infinite' 
% medium, with the initial condition that all mass at $t = 0$ is 
% homogeneously confined inside a sphere of radius $a$. This is sometimes
% called a 'spherical initial condition'.
%
% see J. Crank (1975) "The Mathematics of Diffusion", 2nd Ed., 
%      Clarendon Press (Oxford), pages 29-30 
%      Equation 3.8, Figure 3.1
%
% The transient diffusion equation reads
%
% $$\alpha\frac{\partial c}{\partial t}+\nabla.\left(-D\nabla c\right)=0,$$
%
% where $c$ is the independent variable (concentration, temperature, etc),
% $D$ is the diffusion coefficient, and $\alpha$ is a constant (1, here).
%
%
clc; clear;
more off;
run('../../../FVToolStartUp.m')

%% Define the domain and create a mesh structure
% Here we work in a 1D spherical coordinate system (r coordinate)
L = 10.0;  % domain length
Nx = 2000; % number of cells
m = createMeshSpherical1D(Nx, L);

%% Create the boundary condition structure
BC = createBC(m); % all Neumann boundary condition structure

%% define the transfer coeffs
D = createCellVariable(m, 1.0);
alfa = createCellVariable(m, 1.0);

%% define initial condition
c_init = 0;
c_old = createCellVariable(m, c_init, BC); % initial values
r = c_old.domain.cellcenters.x;
c_old.value(r<1.0) = 1.0;

%% calculate volumes of FV cellslices
%  We use this for demonstrating mass conservation
cellA = m.facecenters.x(1:end-1);
cellB = m.facecenters.x(2:end);
cellvol = 4/3 .* pi .* (cellB.^3 - cellA.^3);
cellsum = sum(cellvol)

c = c_old; % assign the old value of the cells to the current values

t = 0.0; % master time
deltat = 0.0625/20; % time step

% output total mass in the system
m_tot = sum(c.value(2:end-1) .* cellvol);
t,m_tot

%% loop for "time-stepping" the solution
% It outputs the spatial profile C(r) after
% 20, 80 and 320 time-steps
% This corresponds to t=0.0625, t=0.25 and t=1, respectively.
ti = 0
for s=[20,60,240]
  for n=1:s
      [M_trans, RHS_trans] = transientTerm(c, deltat, alfa);
      Dave = harmonicMean(D);
      Mdiff = diffusionTerm(Dave);
      [Mbc, RHSbc] = boundaryCondition(BC);
      M = M_trans-Mdiff+Mbc;
      RHS = RHS_trans+RHSbc;
      c = solvePDE(m,M, RHS);
      t += deltat;
      c_old = c;
  endfor
  m_tot = sum(c.value(2:end-1) .* cellvol);
  n,t,m_tot
  % The following writes the result to a file
  %  adapted from visualizeCells for a spherical 1D mesh (MeshSpherical1D)
  x = [c.domain.facecenters.x(1); c.domain.cellcenters.x; c.domain.facecenters.x(end)];
  cval = [0.5*(c.value(1)+c.value(2)); c.value(2:end-1); 0.5*(c.value(end-1)+c.value(end))];
  ti += s;
  filename = ["diffusion1Dspherical_FVTool_tstep",num2str(ti),".mat"]
  save('-6',filename,'x','cval');
endfor

