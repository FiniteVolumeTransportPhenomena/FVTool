function faceGrad = gradientTermSpherical3D(phi)
% this function calculates the gradient of a variable on a spherical
% (r, theta, phi) 3D mesh, returning a face variable. It uses the ghost
% cell values stored in phi.value.
%
% SYNOPSIS:
%   faceGrad = gradientTermSpherical3D(phi)
%
% SEE ALSO: gradientTermCylindrical3D

% check the size of the variable and the mesh dimension
Nr = phi.domain.dims(1);
Ntheta = phi.domain.dims(2);
Nphi = phi.domain.dims(3);
DR = repmat(phi.domain.cellsize.x, 1, Ntheta, Nphi);
DTHETA = repmat(phi.domain.cellsize.y', Nr, 1, Nphi);
DPHI = zeros(1,1,Nphi+2);
DPHI(1,1,:) = phi.domain.cellsize.z;
DPHI = repmat(DPHI, Nr, Ntheta, 1);
dr = 0.5*(DR(1:end-1,:,:)+DR(2:end,:,:));
dtheta = 0.5*(DTHETA(:,1:end-1,:)+DTHETA(:,2:end,:));
dphi = 0.5*(DPHI(:,:,1:end-1)+DPHI(:,:,2:end));
% metric arrays evaluated on the theta-faces and phi-faces respectively
rp_theta = repmat(phi.domain.cellcenters.x, 1, Ntheta+1, Nphi);
rp_phi = repmat(phi.domain.cellcenters.x, 1, Ntheta, Nphi+1);
thetap_phi = repmat(phi.domain.cellcenters.y', Nr, 1, Nphi+1);

xvalue = (phi.value(2:Nr+2,2:Ntheta+1,2:Nphi+1)-phi.value(1:Nr+1,2:Ntheta+1,2:Nphi+1))./dr;
yvalue = (phi.value(2:Nr+1,2:Ntheta+2,2:Nphi+1)-phi.value(2:Nr+1,1:Ntheta+1,2:Nphi+1))./(dtheta.*rp_theta);
zvalue = (phi.value(2:Nr+1,2:Ntheta+1,2:Nphi+2)-phi.value(2:Nr+1,2:Ntheta+1,1:Nphi+1))./(dphi.*rp_phi.*sin(thetap_phi));

faceGrad=FaceVariable(phi.domain, xvalue, yvalue, zvalue);
end
