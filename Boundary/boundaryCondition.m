function [BCMatrix, BCRHS] = boundaryCondition(BC)
% creates the matrix of coefficient and RHS vector
%
% SYNOPSIS:
%   [BCMatrix, BCRHS] = boundaryCondition(BC)
%
% PARAMETERS:
%   BC             - BoundaryCondition object created by createBC
%
% RETURNS:
%   BCMatrix  - an square sparse matrix
%   BCRHS     - a column vector of values
%
% EXAMPLE:
%   L = 1.0; % length of a 1D domain
%   Nx = 5; % Number of grids in x direction
%   m = createMesh1D(Nx, L);
%   BC = createBC(m); % all Neumann boundaries
%   [Mbc, RHSbc] = boundaryCondition(BC);
%   spy(Mbc); % see inside the boundary matrix of coeffiecients
% SEE ALSO:
%     createBC, createMesh1D, createMesh2D, createMesh3D,
%     createMeshCylindrical1D, createMeshCylindrical2D,
%     createMeshRadial2D, createMeshCylindrical3D,
%     cellBoundary, combineBC, createCellVariable

switch geometryTag(BC.domain)
    case {'1D', 'Cylindrical1D', 'Spherical1D'}
        [BCMatrix, BCRHS] = boundaryCondition1D(BC);
    case {'2D', 'Cylindrical2D'}
        [BCMatrix, BCRHS] = boundaryCondition2D(BC);
    case 'Radial2D'
        [BCMatrix, BCRHS] = boundaryConditionRadial2D(BC);
    case '3D'
        [BCMatrix, BCRHS] = boundaryCondition3D(BC);
    case 'Cylindrical3D'
        [BCMatrix, BCRHS] = boundaryConditionCylindrical3D(BC);
    case 'Spherical3D'
        [BCMatrix, BCRHS] = boundaryConditionSpherical3D(BC);
    otherwise
        error('FVTool:unsupportedGeometry', ...
            'boundaryCondition: no implementation for %s', geometryTag(BC.domain));
end
