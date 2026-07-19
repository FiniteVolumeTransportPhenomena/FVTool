function [M, Mx, My, Mz] = diffusionTerm(D)
% This function uses the central difference scheme to discretize a 2D
% diffusion term in the form \grad . (D \grad \phi) where u is a face vactor
% It also returns the x and y parts of the matrix of coefficient.
%
% SYNOPSIS:
%   [M, Mx, My, Mz] = diffusionTerm(D)
%
% PARAMETERS:
%   D   - diffusion coefficient, FaceVariable
%
% RETURNS:
%  M   - sparse matrix representing the discretized diffusion term
% optionally returns the x and y and z components of the discretized diffusion
% term
% EXAMPLE:
%  m = createMesh2D(3,4,1,1);
%  D = createCellVariable(m,1);
%  D_face = harmonicMean(D);
%  Mdiff = diffusionTerm(D_face);
% 
% SEE ALSO: convectionTerm, createFaceVariable, createCellVariable
%


switch geometryTag(D.domain)
    case '1D'
        M = diffusionTerm1D(D);
    case 'Cylindrical1D'
        M = diffusionTermCylindrical1D(D);
    case 'Spherical1D'
        M = diffusionTermSpherical1D(D);
    case '2D'
        [M, Mx, My] = diffusionTerm2D(D);
    case 'Cylindrical2D'
        [M, Mx, My] = diffusionTermCylindrical2D(D);
    case 'Radial2D'
        [M, Mx, My] = diffusionTermRadial2D(D);
    case '3D'
        [M, Mx, My, Mz] = diffusionTerm3D(D);
    case 'Cylindrical3D'
        [M, Mx, My, Mz] = diffusionTermCylindrical3D(D);
    case 'Spherical3D'
        [M, Mx, My, Mz] = diffusionTermSpherical3D(D);
    otherwise
        error('FVTool:unsupportedGeometry', ...
            'diffusionTerm: no implementation for %s', geometryTag(D.domain));
end
