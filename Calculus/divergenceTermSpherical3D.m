function [RHSdiv, RHSdivx, RHSdivy, RHSdivz] = divergenceTermSpherical3D(F)
% This function calculates the divergence of a field using its face
% average flux vector F on a spherical (r, theta, phi) 3D mesh.
%
% SYNOPSIS:
%   [RHSdiv, RHSdivx, RHSdivy, RHSdivz] = divergenceTermSpherical3D(F)
%
% SEE ALSO: divergenceTermCylindrical3D, divergenceTermSpherical1D

% extract data from the mesh structure
Nx = F.domain.dims(1);
Ntheta = F.domain.dims(2);
Nphi = F.domain.dims(3);
G=reshape((1:(Nx+2)*(Ntheta+2)*(Nphi+2)), Nx+2, Ntheta+2, Nphi+2);
dx = repmat(F.domain.cellsize.x(2:end-1), 1, Ntheta, Nphi);
dy = repmat(F.domain.cellsize.y(2:end-1)', Nx, 1, Nphi);
DZ = zeros(1,1,Nphi+2);
DZ(1,1,:) = F.domain.cellsize.z;
dz=repmat(DZ(1,1,2:end-1), Nx, Ntheta, 1);
rp = repmat(F.domain.cellcenters.x, 1, Ntheta, Nphi);
rf = repmat(F.domain.facecenters.x, 1, Ntheta, Nphi);
thetap = repmat(F.domain.cellcenters.y', Nx, 1, Nphi);
thetaf = repmat(F.domain.facecenters.y', Nx, 1, Nphi);

% define the vector of cell index
row_index = reshape(G(2:Nx+1,2:Ntheta+1, 2:Nphi+1),Nx*Ntheta*Nphi,1); % main diagonal (only internal cells)

% reassign the flux vectors for code readability
Fx = F.xvalue;
Fy = F.yvalue;
Fz = F.zvalue;

% compute the divergence with the spherical metric factors
div_x = (rf(2:Nx+1,:,:).^2.*Fx(2:Nx+1,:,:) - rf(1:Nx,:,:).^2.*Fx(1:Nx,:,:))./(dx.*rp.^2);
div_y = (sin(thetaf(:,2:Ntheta+1,:)).*Fy(:,2:Ntheta+1,:) - sin(thetaf(:,1:Ntheta,:)).*Fy(:,1:Ntheta,:))./(dy.*rp.*sin(thetap));
div_z = (Fz(:,:,2:Nphi+1) - Fz(:,:,1:Nphi))./(dz.*rp.*sin(thetap));

% define the RHS Vector
RHSdiv = zeros((Nx+2)*(Ntheta+2)*(Nphi+2),1);
RHSdivx = zeros((Nx+2)*(Ntheta+2)*(Nphi+2),1);
RHSdivy = zeros((Nx+2)*(Ntheta+2)*(Nphi+2),1);
RHSdivz = zeros((Nx+2)*(Ntheta+2)*(Nphi+2),1);

% assign the values of the RHS vector
RHSdiv(row_index) = reshape(div_x+div_y+div_z,Nx*Ntheta*Nphi,1);
RHSdivx(row_index) = reshape(div_x,Nx*Ntheta*Nphi,1);
RHSdivy(row_index) = reshape(div_y,Nx*Ntheta*Nphi,1);
RHSdivz(row_index) = reshape(div_z,Nx*Ntheta*Nphi,1);
end
