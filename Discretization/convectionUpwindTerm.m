function [M, Mx, My, Mz] = convectionUpwindTerm(u, varargin)
% This function uses the upwind scheme to discretize a 2D
% convection term in the form \grad (u \phi) where u is a face vactor
% It also returns the x and y parts of the matrix of coefficient.
%
% SYNOPSIS:
%   [M, Mx, My, Mz] = convectionUpwindTerm(u)
%
% PARAMETERS:
%   u   - velovity vector, face variable
%
% RETURNS:
%
%
% EXAMPLE:
%
% SEE ALSO:
%

if nargin>1
    switch geometryTag(u.domain)
        case '1D'
            M = convectionUpwindTerm1D(u, varargin{1});
        case 'Cylindrical1D'
            M = convectionUpwindTermCylindrical1D(u, varargin{1});
        case 'Spherical1D'
            M = convectionUpwindTermSpherical1D(u, varargin{1});
        case '2D'
            [M, Mx, My] = convectionUpwindTerm2D(u, varargin{1});
        case 'Cylindrical2D'
            [M, Mx, My] = convectionUpwindTermCylindrical2D(u, varargin{1});
        case 'Radial2D'
            [M, Mx, My] = convectionUpwindTermRadial2D(u, varargin{1});
        case '3D'
            [M, Mx, My, Mz] = convectionUpwindTerm3D(u, varargin{1});
        case 'Cylindrical3D'
            [M, Mx, My, Mz] = convectionUpwindTermCylindrical3D(u, varargin{1});
        case 'Spherical3D'
            [M, Mx, My, Mz] = convectionUpwindTermSpherical3D(u, varargin{1});
        otherwise
            error('FVTool:unsupportedGeometry', ...
                'convectionUpwindTerm: no implementation for %s', geometryTag(u.domain));
    end
else
    switch geometryTag(u.domain)
        case '1D'
            M = convectionUpwindTerm1D(u);
        case 'Cylindrical1D'
            M = convectionUpwindTermCylindrical1D(u);
        case 'Spherical1D'
            M = convectionUpwindTermSpherical1D(u);
        case '2D'
            [M, Mx, My] = convectionUpwindTerm2D(u);
        case 'Cylindrical2D'
            [M, Mx, My] = convectionUpwindTermCylindrical2D(u);
        case 'Radial2D'
            [M, Mx, My] = convectionUpwindTermRadial2D(u);
        case '3D'
            [M, Mx, My, Mz] = convectionUpwindTerm3D(u);
        case 'Cylindrical3D'
            [M, Mx, My, Mz] = convectionUpwindTermCylindrical3D(u);
        case 'Spherical3D'
            [M, Mx, My, Mz] = convectionUpwindTermSpherical3D(u);
        otherwise
            error('FVTool:unsupportedGeometry', ...
                'convectionUpwindTerm: no implementation for %s', geometryTag(u.domain));
    end
end