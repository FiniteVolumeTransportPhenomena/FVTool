function [RHSdiv, RHSdivx, RHSdivy, RHSdivz] = divergenceTerm(F)
% This function calculates the divergence of a field using its face
% average value and the vector F, which is a face vector
%
% SYNOPSIS:
%   [RHSdiv, RHSdivx, RHSdivy, RHSdivz] = divergenceTerm(F)
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

RHSdivz=[];

switch geometryTag(F.domain)
    case '1D'
        RHSdiv = divergenceTerm1D(F);
    case 'Cylindrical1D'
        RHSdiv = divergenceTermCylindrical1D(F);
    case 'Spherical1D'
        RHSdiv = divergenceTermSpherical1D(F);
    case '2D'
        [RHSdiv, RHSdivx, RHSdivy] = divergenceTerm2D(F);
    case 'Cylindrical2D'
        [RHSdiv, RHSdivx, RHSdivy] = divergenceTermCylindrical2D(F);
    case 'Radial2D'
        [RHSdiv, RHSdivx, RHSdivy] = divergenceTermRadial2D(F);
    case '3D'
        [RHSdiv, RHSdivx, RHSdivy, RHSdivz] = divergenceTerm3D(F);
    case 'Cylindrical3D'
        [RHSdiv, RHSdivx, RHSdivy, RHSdivz] = divergenceTermCylindrical3D(F);
    case 'Spherical3D'
        [RHSdiv, RHSdivx, RHSdivy, RHSdivz] = divergenceTermSpherical3D(F);
    otherwise
        error('FVTool:unsupportedGeometry', ...
            'divergenceTerm: no implementation for %s', geometryTag(F.domain));
end
