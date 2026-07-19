function phiBC = cellBoundary(phi, BC)
% This function calculates the value of the boundary cells and add them
% to the variable phi. The function is used to add ghost cells to a matrix
% that specifies the cell values over a domain
%
% SYNOPSIS:
%   phiBC = cellBoundary(phi, BC)
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

% extract data from the mesh structure
switch geometryTag(BC.domain)
    case {'1D', 'Cylindrical1D', 'Spherical1D'}
        phiBC = cellBoundary1D(phi, BC);
    case {'2D', 'Cylindrical2D'}
        phiBC = cellBoundary2D(phi, BC);
    case 'Radial2D'
        phiBC = cellBoundaryRadial2D(phi, BC);
    case '3D'
        phiBC = cellBoundary3D(phi, BC);
    case 'Cylindrical3D'
        phiBC = cellBoundaryCylindrical3D(phi, BC);
    case 'Spherical3D'
        phiBC = cellBoundarySpherical3D(phi, BC);
    otherwise
        error('FVTool:unsupportedGeometry', ...
            'cellBoundary: no implementation for %s', geometryTag(BC.domain));
end
