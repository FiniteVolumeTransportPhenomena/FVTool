function BC = createBC(meshvar)
% function BC = createBC(meshvar)
% Creates a boundary condition structure from a mesh structure
% The boundary conditions on all boundaries are Neumann by default
% and can be altered later bu the user; see examples
%
% SYNOPSIS:
%	BC = createBC(meshvar)
%
% PARAMETERS:
%
%
% RETURNS:
%
%
% EXAMPLE:
%	L = 1.0; % length of a 1D domain
%   Nx = 5; % Number of grids in x direction
%   m = createMesh1D(Nx, L);
%   BC = createBC(m); % all Neumann boundaries
%
% SEE ALSO:
%

switch meshvar.dimension
    case 1
        BC = createBC1D(meshvar);
    case 2
        BC = createBC2D(meshvar);
    case 3
        BC = createBC3D(meshvar);
    otherwise
        error('FVTool:unsupportedGeometry', ...
            'createBC: unsupported dimension %g', meshvar.dimension);
end
