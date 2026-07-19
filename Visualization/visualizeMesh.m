function visualizeMesh(MS)
%VISUALIZECELLS plots the values of cell variable phi.value
%
% SYNOPSIS:
%   visualizeCells(MS)
%
% PARAMETERS:
%  MS: Mesh Structure
%
% RETURNS:
%    Nothing
%
% EXAMPLE:
%
% SEE ALSO:
%

% Written by Ali A. Eftekhari
% See the license file

switch geometryTag(MS)
    case '1D'
        plot(MS.facecenters.x, zeros(size(MS.facecenters.x)), '-+');
        title('1D Cartesian Grid')
    case 'Cylindrical1D'
        L = MS.facecenters.x(end);
        [TH,R] = meshgrid(linspace(-pi/16, pi/16, 3), MS.facecenters.x);
        [X,Y] = pol2cart(TH,R);
        h = polar([0 2*pi], [0 L]);
        delete(h);
        hold on
        pcolor(X,Y,zeros(size(X)));
        title('1D Radial Grid')
        hold off
    case {'2D', 'Cylindrical2D'}
        phi.value(:,1) = 0.5*(phi.value(:,1)+phi.value(:,2));
        phi.value(1,:) = 0.5*(phi.value(1,:)+phi.value(2,:));
        phi.value(:,end) = 0.5*(phi.value(:,end)+phi.value(:,end-1));
        phi.value(end,:) = 0.5*(phi.value(end,:)+phi.value(end-1,:));
        phi.value(1,1) = phi.value(1,2); phi.value(1,end) = phi.value(1,end-1);
        phi.value(end,1) = phi.value(end,2); phi.value(end,end) = phi.value(end,end-1);
        visualizeCells2D(phi);
    case 'Radial2D'
        [TH,R] = meshgrid(MS.facecenters.y, MS.facecenters.x);
        [X,Y] = pol2cart(TH,R);
        h = polar([0 2*pi], [0 L]);
        delete(h);
        hold on
        pcolor(X,Y,phi.value)
        colorbar
        hold off
    case '3D'
        phi.value = phi.value(2:end-1,2:end-1,2:end-1);
        visualizeCells3D(phi);
    case 'Cylindrical3D'
        phi.value = phi.value(2:end-1,2:end-1,2:end-1);
        visualizeCellsCylindrical3D(phi);
    otherwise
        error('FVTool:unsupportedGeometry', ...
            'visualizeMesh: no implementation for %s', geometryTag(MS));
end

end
