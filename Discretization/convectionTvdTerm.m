function [M, RHS, Mx, My, Mz, RHSx, RHSy, RHSz] = convectionTvdTerm(u, phi, FL)
% This function uses the TVD scheme to discretize a
% convection term in the form $\grad (u \phi)$ where u is a face vactor
% It also returns the x, y, x parts of the matrix of coefficient.
%
% SYNOPSIS:
%   [M, RHS, Mx, My, Mz, RHSx, RHSy, RHSz] = convectionTvdTerm(u, phi, FL)
%
% PARAMETERS:
%   u  - velocity vector, FaceVariable
%   phi  - value of phi from the previous time step or iteration, CellVariable
%   FL  - Flux Limiter function
%
% RETURNS:
%
%
% EXAMPLE:
%
% SEE ALSO:
%

Mz=[];
switch geometryTag(u.domain)
    case '1D'
        [M, RHS] = convectionTvdTerm1D(u, phi, FL);
    case 'Cylindrical1D'
        [M, RHS] = convectionTvdTermCylindrical1D(u, phi, FL);
    case 'Spherical1D'
        [M, RHS] = convectionTvdTermSpherical1D(u, phi, FL);
    case '2D'
        [M, RHS, Mx, My, RHSx, RHSy] = convectionTvdTerm2D(u, phi, FL);
    case 'Cylindrical2D'
        [M, RHS, Mx, My, RHSx, RHSy] = convectionTvdTermCylindrical2D(u, phi, FL);
    case 'Radial2D'
        [M, RHS, Mx, My, RHSx, RHSy] = ...
            convectionTvdTermRadial2D(u, phi, FL);
    case '3D'
        [M, RHS, Mx, My, Mz, RHSx, RHSy, RHSz] = convectionTvdTerm3D(u, phi, FL);
    case 'Cylindrical3D'
        [M, RHS, Mx, My, Mz, RHSx, RHSy, RHSz] = ...
            convectionTvdTermCylindrical3D(u, phi, FL);
    case 'Spherical3D'
        [M, RHS, Mx, My, Mz, RHSx, RHSy, RHSz] = ...
            convectionTvdTermSpherical3D(u, phi, FL);
    otherwise
        error('FVTool:unsupportedGeometry', ...
            'convectionTvdTerm: no implementation for %s', geometryTag(u.domain));
end
