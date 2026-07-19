function visualizeCellVectors(phi_cell)
%VISUALIZECELLS plots the values of cell variable phi.value
%
% SYNOPSIS:
%  visualizeCellVectors(phi_cell)
%
% PARAMETERS:
%  phi_cell: a CellVector variable
%
% RETURNS:
%  None
%
% EXAMPLE:
%
% SEE ALSO:
%

% Written by Ali A. Eftekhari
% See the license file

switch geometryTag(phi_cell.domain)
    case {'1D', 'Cylindrical1D', 'Spherical1D'}
        warning('No vector visualization for a 1D domain.');
    case {'2D', 'Cylindrical2D'}
        visualizeCellVectors2D(phi_cell);
    case 'Radial2D'
        visualizeCellVectorsRadial2D(phi_cell);
    case '3D'
        visualizeCellVectors3D(phi_cell);
    case 'Cylindrical3D'
        visualizeCellVectorsCylindrical3D(phi_cell);
    otherwise
        error('FVTool:unsupportedGeometry', ...
            'visualizeCellVectors: no implementation for %s', geometryTag(phi_cell.domain));
end

end
