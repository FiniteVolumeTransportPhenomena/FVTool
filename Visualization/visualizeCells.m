function visualizeCells(phi)
%VISUALIZECELLS plots the values of cell variable phi.value
%
% SYNOPSIS:
%
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

% Written by Ali A. Eftekhari
% See the license file

switch geometryTag(phi.domain)
    case {'1D', 'Cylindrical1D', 'Spherical1D'}
        phi.value = [0.5*(phi.value(1)+phi.value(2)); phi.value(2:end-1); 0.5*(phi.value(end-1)+phi.value(end))];
        visualizeCells1D(phi);
    case {'2D', 'Cylindrical2D', 'Radial2D'}
        phi.value(:,1) = 0.5*(phi.value(:,1)+phi.value(:,2));
        phi.value(1,:) = 0.5*(phi.value(1,:)+phi.value(2,:));
        phi.value(:,end) = 0.5*(phi.value(:,end)+phi.value(:,end-1));
        phi.value(end,:) = 0.5*(phi.value(end,:)+phi.value(end-1,:));
        phi.value(1,1) = phi.value(1,2); phi.value(1,end) = phi.value(1,end-1);
        phi.value(end,1) = phi.value(end,2); phi.value(end,end) = phi.value(end,end-1);
        if strcmp(geometryTag(phi.domain), 'Radial2D')
            visualizeCellsRadial2D(phi);
        else
            visualizeCells2D(phi);
        end
    case '3D'
        phi.value = phi.value(2:end-1,2:end-1,2:end-1);
        visualizeCells3D(phi);
    case 'Cylindrical3D'
        phi.value = phi.value(2:end-1,2:end-1,2:end-1);
        visualizeCellsCylindrical3D(phi);
    case 'Spherical3D'
        phi.value = phi.value(2:end-1,2:end-1,2:end-1);
        visualizeCellsSpherical3D(phi);
    otherwise
        error('FVTool:unsupportedGeometry', ...
            'visualizeCells: no implementation for %s', geometryTag(phi.domain));
end

end
