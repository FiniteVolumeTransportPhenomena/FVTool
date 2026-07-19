function visualizeCellsSpherical3D(phi)
%VISUALIZECELLSSPHERICAL3D plots the values of a cell variable phi defined
% on a spherical (r, theta, phi) 3D mesh.
%
% SYNOPSIS:
%   visualizeCellsSpherical3D(phi)
%
% SEE ALSO: visualizeCellsCylindrical3D

% Written by Ali A. Eftekhari
% See the license file

d = phi.domain.dims;
[TH,R,PH]=meshgrid(phi.domain.cellcenters.y, ...
    phi.domain.cellcenters.x, phi.domain.cellcenters.z);
% spherical (r, theta measured from +z axis, phi azimuth) to cartesian
X = R.*sin(TH).*cos(PH);
Y = R.*sin(TH).*sin(PH);
Z = R.*cos(TH);
hold on
% surf three spherical shells (constant r)
surf(squeeze(X(floor(2*d(1)/3),:,:)), squeeze(Y(floor(2*d(1)/3),:,:)), ...
    squeeze(Z(floor(2*d(1)/3),:,:)), squeeze(phi.value(floor(2*d(1)/3),:,:)));
surf(squeeze(X(floor(d(1)/3),:,:)), squeeze(Y(floor(d(1)/3),:,:)), ...
    squeeze(Z(floor(d(1)/3),:,:)), squeeze(phi.value(floor(d(1)/3),:,:)));
surf(squeeze(X(1,:,:)), squeeze(Y(1,:,:)), ...
    squeeze(Z(1,:,:)), squeeze(phi.value(1,:,:)));

% surf constant-phi meridional sections
surf(squeeze(X(:,:,1)), squeeze(Y(:,:,1)), ...
    squeeze(Z(:,:,1)), squeeze(phi.value(:,:,1)));
surf(squeeze(X(:,:,floor(d(3)/2))), squeeze(Y(:,:,floor(d(3)/2))), ...
    squeeze(Z(:,:,floor(d(3)/2))), squeeze(phi.value(:,:,floor(d(3)/2))));
surf(squeeze(X(:,:,end)), squeeze(Y(:,:,end)), ...
    squeeze(Z(:,:,end)), squeeze(phi.value(:,:,end)));

% surf constant-theta conical sections
surf(squeeze(X(:,floor(d(2)/3),:)), squeeze(Y(:,floor(d(2)/3),:)), ...
    squeeze(Z(:,floor(d(2)/3),:)), squeeze(phi.value(:,floor(d(2)/3),:)));
surf(squeeze(X(:,floor(2*d(2)/3),:)), squeeze(Y(:,floor(2*d(2)/3),:)), ...
    squeeze(Z(:,floor(2*d(2)/3),:)), squeeze(phi.value(:,floor(2*d(2)/3),:)));

axis equal tight
view([60 25]);
colorbar;
hold off
end
