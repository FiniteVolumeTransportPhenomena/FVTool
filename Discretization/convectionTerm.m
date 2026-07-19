function [M, Mx, My, Mz] = convectionTerm(u)
% This function uses the central difference scheme to discretize a 2D
% convection term in the form \grad (u \phi) where u is a face vactor
% It also returns the x and y parts of the matrix of coefficient.
%
% SYNOPSIS:
%   [M, Mx, My, Mz] = convectionTerm(u)
%
% PARAMETERS:
%   u   - FaceVariable  
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
    	M = convectionTerm1D(u);
    case 'Cylindrical1D'
        M = convectionTermCylindrical1D(u);
    case 'Spherical1D'
        M = convectionTermSpherical1D(u);
    case '2D'
        [M, Mx, My] = convectionTerm2D(u);
    case 'Cylindrical2D'
        [M, Mx, My] = convectionTermCylindrical2D(u);
    case 'Radial2D'
        [M, Mx, My] = convectionTermRadial2D(u);
    case '3D'
        [M, Mx, My, Mz] = convectionTerm3D(u);
    case 'Cylindrical3D'
        [M, Mx, My, Mz] = convectionTermCylindrical3D(u);
    case 'Spherical3D'
        [M, Mx, My, Mz] = convectionTermSpherical3D(u);
    otherwise
        error('FVTool:unsupportedGeometry', ...
            'convectionTerm: no implementation for %s', geometryTag(u.domain));
end
