function [RHS, RHSx, RHSy, RHSz] = convectionTvdRHS(u, phi, FL)
% This function uses the TVD scheme to discretize a
% convection term in the form $\grad (u \phi)$ where u is a face vactor
% It also returns the x, y, x parts of the matrix of coefficient.
%
% SYNOPSIS:
%   [RHS, RHSx, RHSy, RHSz] = convectionTvdRHS(u, phi, FL)
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

switch geometryTag(u.domain)
    case '1D'
        RHS = convectionTvdRHS1D(u, phi, FL);
    case 'Cylindrical1D'
        RHS = convectionTvdRHSCylindrical1D(u, phi, FL);
    case 'Spherical1D'
        RHS = convectionTvdRHSSpherical1D(u, phi, FL);
    case '2D'
        [RHS, RHSx, RHSy] = convectionTvdRHS2D(u, phi, FL);
    case 'Cylindrical2D'
        [RHS, RHSx, RHSy] = convectionTvdRHSCylindrical2D(u, phi, FL);
    case 'Radial2D'
        [RHS, RHSx, RHSy] = ...
            convectionTvdRHSRadial2D(u, phi, FL);
    case '3D'
        [RHS, RHSx, RHSy, RHSz] = convectionTvdRHS3D(u, phi, FL);
    case 'Cylindrical3D'
        [RHS, RHSx, RHSy, RHSz] = ...
            convectionTvdRHSCylindrical3D(u, phi, FL);
    case 'Spherical3D'
        [RHS, RHSx, RHSy, RHSz] = ...
            convectionTvdRHSSpherical3D(u, phi, FL);
    otherwise
        error('FVTool:unsupportedGeometry', ...
            'convectionTvdRHS: no implementation for %s', geometryTag(u.domain));
end
