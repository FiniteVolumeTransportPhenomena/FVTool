function phi = solveExplicitPDE(phi_old, dt, RHS, BC, varargin)
% phi = solveExplicitPDE(phi_old, dt, RHS, BC, varargin)
% SolveExplicitPDE solves for the new value of variable \phi in an explicit
% discretization scheme: \phi_new = \phi_old + dt/ * RHS
% the real equation is d\phi/dt=RHS
% The code calculates the new values for the internal cells. Then it uses
% the boundary condition to calculate the values for the ghost cells

% SYNOPSIS:
%
%
% PARAMETERS:
%
%
% RETURNS:
%
%
% EXAMPLE:
%
% SEE ALSO:
%

% Written by Ali A. Eftekhari
% See the license file

% this is not the most beautiful implementation, but it works.
if nargin==5
    x = phi_old.value(:)+dt*RHS./varargin{1}.value(:);
else
    x = phi_old.value(:)+dt*RHS;
end

d = phi_old.domain.dimension;
N = phi_old.domain.dims;

if (d>=2)
    phi_val = reshape(x, N+2);
else
    phi_val = reshape(x, [N(1)+2 1]);
end

if (d == 1)
	phi_val = phi_val(2:N(1)+1);
elseif (d == 2)
	phi_val = phi_val(2:N(1)+1, 2:N(2)+1);
elseif (d == 3)
    phi_val = phi_val(2:N(1)+1, 2:N(2)+1, 2:N(3)+1);
end

phi= createCellVariable(phi_old.domain, phi_val, BC);
end
