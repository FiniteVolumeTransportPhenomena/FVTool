function faceGrad = gradientTerm(phi)
% this function calculates the gradient of a variable in x direction
% it checks for the availability of the ghost variables and use them, otherwise
% estimate them, assuming a zero gradient on the boundaries
%
% SYNOPSIS:
%   faceGrad = gradientTerm(phi)
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

% Copyright (c) 2012-2016 Ali Akbar Eftekhari
% See the license file

switch geometryTag(phi.domain)
    case {'1D', 'Cylindrical1D', 'Spherical1D'}
        faceGrad = gradientTerm1D(phi);
    case {'2D', 'Cylindrical2D'}
        faceGrad = gradientTerm2D(phi);
    case 'Radial2D'
        faceGrad = gradientTermRadial2D(phi);
    case '3D'
        faceGrad = gradientTerm3D(phi);
    case 'Cylindrical3D'
        faceGrad = gradientTermCylindrical3D(phi);
    case 'Spherical3D'
        faceGrad = gradientTermSpherical3D(phi);
    otherwise
        error('FVTool:unsupportedGeometry', ...
            'gradientTerm: no implementation for %s', geometryTag(phi.domain));
end
