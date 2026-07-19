function [M, Mx, My, Mz] = diffusionTermSpherical3D(D)
% This function uses the central difference scheme to discretize a 3D
% diffusion term in the form \grad . (D \grad \phi) on a spherical
% (r, theta, phi) grid, where D is a face variable.
% It also returns the r (x), theta (y) and phi (z) parts of the matrix
% of coefficients.
%
% SYNOPSIS:
%   [M, Mx, My, Mz] = diffusionTermSpherical3D(D)
%
% PARAMETERS:
%   D   - diffusion coefficient, FaceVariable on a MeshSpherical3D
%
% RETURNS:
%  M   - sparse matrix representing the discretized diffusion term
%  optionally the r, theta and phi components of the term.
%
% SEE ALSO: diffusionTermCylindrical3D, diffusionTermSpherical1D

% extract data from the mesh structure
Nr = D.domain.dims(1);
Ntheta = D.domain.dims(2);
Nphi = D.domain.dims(3);
G=reshape(1:(Nr+2)*(Ntheta+2)*(Nphi+2), Nr+2, Ntheta+2, Nphi+2);
DR = repmat(D.domain.cellsize.x, 1, Ntheta, Nphi);
DTHETA = repmat(D.domain.cellsize.y', Nr, 1, Nphi);
DPHI = ones(1,1,Nphi+2);
DPHI(1,1,:) = D.domain.cellsize.z;
DPHI = repmat(DPHI, Nr, Ntheta, 1);
dr = 0.5*(DR(1:end-1,:,:)+DR(2:end,:,:));
dtheta = 0.5*(DTHETA(:,1:end-1,:)+DTHETA(:,2:end,:));
dphi = 0.5*(DPHI(:,:,1:end-1)+DPHI(:,:,2:end));
rp = repmat(D.domain.cellcenters.x, 1, Ntheta, Nphi);
rf = repmat(D.domain.facecenters.x, 1, Ntheta, Nphi);
thetap = repmat(D.domain.cellcenters.y', Nr, 1, Nphi);
thetaf = repmat(D.domain.facecenters.y', Nr, 1, Nphi);

% define the vectors to stores the sparse matrix data
iix = zeros(3*(Nr+2)*(Ntheta+2)*(Nphi+2),1);
jjx = zeros(3*(Nr+2)*(Ntheta+2)*(Nphi+2),1);
sx = zeros(3*(Nr+2)*(Ntheta+2)*(Nphi+2),1);
iiy = zeros(3*(Nr+2)*(Ntheta+2)*(Nphi+2),1);
jjy = zeros(3*(Nr+2)*(Ntheta+2)*(Nphi+2),1);
sy = zeros(3*(Nr+2)*(Ntheta+2)*(Nphi+2),1);
iiz = zeros(3*(Nr+2)*(Ntheta+2)*(Nphi+2),1);
jjz = zeros(3*(Nr+2)*(Ntheta+2)*(Nphi+2),1);
sz = zeros(3*(Nr+2)*(Ntheta+2)*(Nphi+2),1);
mNr = Nr*Ntheta*Nphi;	mny = Nr*Ntheta*Nphi;   mnz = Nr*Ntheta*Nphi;

% reassign the east (r+), west (r-), north (theta+), south (theta-),
% front (phi+) and back (phi-) diffusion coefficients following the
% spherical Laplacian metric factors
De = rf(2:Nr+1,:,:).^2.*D.xvalue(2:Nr+1,:,:)./(rp.^2.*dr(2:Nr+1,:,:).*DR(2:Nr+1,:,:));
Dw = rf(1:Nr,:,:).^2.*D.xvalue(1:Nr,:,:)./(rp.^2.*dr(1:Nr,:,:).*DR(2:Nr+1,:,:));
Dn = D.yvalue(:,2:Ntheta+1,:).*sin(thetaf(:,2:Ntheta+1,:))./(rp.^2.*sin(thetap).*dtheta(:,2:Ntheta+1,:).*DTHETA(:,2:Ntheta+1,:));
Ds = D.yvalue(:,1:Ntheta,:).*sin(thetaf(:,1:Ntheta,:))./(rp.^2.*sin(thetap).*dtheta(:,1:Ntheta,:).*DTHETA(:,2:Ntheta+1,:));
Df = D.zvalue(:,:,2:Nphi+1)./(rp.^2.*sin(thetap).^2.*dphi(:,:,2:Nphi+1).*DPHI(:,:,2:Nphi+1));
Db = D.zvalue(:,:,1:Nphi)./(rp.^2.*sin(thetap).^2.*dphi(:,:,1:Nphi).*DPHI(:,:,2:Nphi+1));

% calculate the coefficients for the internal cells
AE = reshape(De,mNr,1);
AW = reshape(Dw,mNr,1);
AN = reshape(Dn,mny,1);
AS = reshape(Ds,mny,1);
AF = reshape(Df,mnz,1);
AB = reshape(Db,mnz,1);
APx = reshape(-(De+Dw),mNr,1);
APy = reshape(-(Dn+Ds),mny,1);
APz = reshape(-(Df+Db),mnz,1);

% build the sparse matrix based on the numbering system
rowx_index = reshape(G(2:Nr+1,2:Ntheta+1,2:Nphi+1),mNr,1); % main diagonal x
iix(1:3*mNr) = repmat(rowx_index,3,1);
rowy_index = reshape(G(2:Nr+1,2:Ntheta+1,2:Nphi+1),mny,1); % main diagonal y
iiy(1:3*mny) = repmat(rowy_index,3,1);
rowz_index = reshape(G(2:Nr+1,2:Ntheta+1,2:Nphi+1),mnz,1); % main diagonal z
iiz(1:3*mnz) = repmat(rowz_index,3,1);
jjx(1:3*mNr) = [reshape(G(1:Nr,2:Ntheta+1,2:Nphi+1),mNr,1); reshape(G(2:Nr+1,2:Ntheta+1,2:Nphi+1),mNr,1); reshape(G(3:Nr+2,2:Ntheta+1,2:Nphi+1),mNr,1)];
jjy(1:3*mny) = [reshape(G(2:Nr+1,1:Ntheta,2:Nphi+1),mny,1); reshape(G(2:Nr+1,2:Ntheta+1,2:Nphi+1),mny,1); reshape(G(2:Nr+1,3:Ntheta+2,2:Nphi+1),mny,1)];
jjz(1:3*mnz) = [reshape(G(2:Nr+1,2:Ntheta+1,1:Nphi),mnz,1); reshape(G(2:Nr+1,2:Ntheta+1,2:Nphi+1),mnz,1); reshape(G(2:Nr+1,2:Ntheta+1,3:Nphi+2),mnz,1)];
sx(1:3*mNr) = [AW; APx; AE];
sy(1:3*mny) = [AS; APy; AN];
sz(1:3*mnz) = [AB; APz; AF];

% build the sparse matrix
kx = 3*mNr;
ky = 3*mny;
kz = 3*mnz;
Mx = sparse(iix(1:kx), jjx(1:kx), sx(1:kx), (Nr+2)*(Ntheta+2)*(Nphi+2), (Nr+2)*(Ntheta+2)*(Nphi+2));
My = sparse(iiy(1:ky), jjy(1:ky), sy(1:ky), (Nr+2)*(Ntheta+2)*(Nphi+2), (Nr+2)*(Ntheta+2)*(Nphi+2));
Mz = sparse(iiz(1:kz), jjz(1:kz), sz(1:kz), (Nr+2)*(Ntheta+2)*(Nphi+2), (Nr+2)*(Ntheta+2)*(Nphi+2));
M = Mx + My + Mz;
end
