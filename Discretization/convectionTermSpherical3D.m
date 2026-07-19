function [M, Mx, My, Mz] = convectionTermSpherical3D(u)
% This function uses the central difference scheme to discretize a 3D
% convection term in the form \grad . (u \phi) on a spherical
% (r, theta, phi) grid, where u is a face vector.
%
% SYNOPSIS:
%   [M, Mx, My, Mz] = convectionTermSpherical3D(u)
%
% SEE ALSO: convectionTermCylindrical3D

% extract data from the mesh structure
Nr = u.domain.dims(1);
Ntheta = u.domain.dims(2);
Nphi = u.domain.dims(3);
G=reshape((1:(Nr+2)*(Ntheta+2)*(Nphi+2)), Nr+2, Ntheta+2, Nphi+2);
DRe = repmat(u.domain.cellsize.x(3:end), 1, Ntheta, Nphi);
DRw = repmat(u.domain.cellsize.x(1:end-2), 1, Ntheta, Nphi);
DRp = repmat(u.domain.cellsize.x(2:end-1), 1, Ntheta, Nphi);
DTHETAn = repmat(u.domain.cellsize.y(3:end)', Nr, 1, Nphi);
DTHETAs = repmat(u.domain.cellsize.y(1:end-2)', Nr, 1, Nphi);
DTHETAp = repmat(u.domain.cellsize.y(2:end-1)', Nr, 1, Nphi);
DPHI = zeros(1,1,Nphi+2);
DPHI(1,1,:) = u.domain.cellsize.z;
DPHIf=repmat(DPHI(1,1,3:end), Nr, Ntheta, 1);
DPHIb=repmat(DPHI(1,1,1:end-2), Nr, Ntheta, 1);
DPHIp=repmat(DPHI(1,1,2:end-1), Nr, Ntheta, 1);
rp = repmat(u.domain.cellcenters.x, 1, Ntheta, Nphi);
rf = repmat(u.domain.facecenters.x, 1, Ntheta, Nphi);
thetap = repmat(u.domain.cellcenters.y', Nr, 1, Nphi);
thetaf = repmat(u.domain.facecenters.y', Nr, 1, Nphi);

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
mnx = Nr*Ntheta*Nphi;	mny = Nr*Ntheta*Nphi;   mnz = Nr*Ntheta*Nphi;

% reassign the east/west (r), north/south (theta), front/back (phi)
% velocity vectors following the spherical metric factors
ue = rf(2:Nr+1,:,:).^2.*u.xvalue(2:Nr+1,:,:)./(rp.^2.*(DRp+DRe));
uw = rf(1:Nr,:,:).^2.*u.xvalue(1:Nr,:,:)./(rp.^2.*(DRp+DRw));
vn = u.yvalue(:,2:Ntheta+1,:).*sin(thetaf(:,2:Ntheta+1,:))./(rp.*sin(thetap).*(DTHETAp+DTHETAn));
vs = u.yvalue(:,1:Ntheta,:).*sin(thetaf(:,1:Ntheta,:))./(rp.*sin(thetap).*(DTHETAp+DTHETAs));
wf = u.zvalue(:,:,2:Nphi+1)./(rp.*sin(thetap).*(DPHIp+DPHIf));
wb = u.zvalue(:,:,1:Nphi)./(rp.*sin(thetap).*(DPHIp+DPHIb));

% calculate the coefficients for the internal cells
AE = reshape(ue,mnx,1);
AW = reshape(-uw,mnx,1);
AN = reshape(vn,mny,1);
AS = reshape(-vs,mny,1);
AF = reshape(wf,mnz,1);
AB = reshape(-wb,mnz,1);
APx = reshape((DRe.*ue-DRw.*uw)./DRp,mnx,1);
APy = reshape((DTHETAn.*vn-DTHETAs.*vs)./DTHETAp,mny,1);
APz = reshape((DPHIf.*wf-DPHIb.*wb)./DPHIp,mnz,1);

% build the sparse matrix based on the numbering system
rowx_index = reshape(G(2:Nr+1,2:Ntheta+1,2:Nphi+1),mnx,1); % main diagonal x
iix(1:3*mnx) = repmat(rowx_index,3,1);
rowy_index = reshape(G(2:Nr+1,2:Ntheta+1,2:Nphi+1),mny,1); % main diagonal y
iiy(1:3*mny) = repmat(rowy_index,3,1);
rowz_index = reshape(G(2:Nr+1,2:Ntheta+1,2:Nphi+1),mnz,1); % main diagonal z
iiz(1:3*mnz) = repmat(rowz_index,3,1);
jjx(1:3*mnx) = [reshape(G(1:Nr,2:Ntheta+1,2:Nphi+1),mnx,1); reshape(G(2:Nr+1,2:Ntheta+1,2:Nphi+1),mnx,1); reshape(G(3:Nr+2,2:Ntheta+1,2:Nphi+1),mnx,1)];
jjy(1:3*mny) = [reshape(G(2:Nr+1,1:Ntheta,2:Nphi+1),mny,1); reshape(G(2:Nr+1,2:Ntheta+1,2:Nphi+1),mny,1); reshape(G(2:Nr+1,3:Ntheta+2,2:Nphi+1),mny,1)];
jjz(1:3*mnz) = [reshape(G(2:Nr+1,2:Ntheta+1,1:Nphi),mnz,1); reshape(G(2:Nr+1,2:Ntheta+1,2:Nphi+1),mnz,1); reshape(G(2:Nr+1,2:Ntheta+1,3:Nphi+2),mnz,1)];
sx(1:3*mnx) = [AW; APx; AE];
sy(1:3*mny) = [AS; APy; AN];
sz(1:3*mnz) = [AB; APz; AF];

% build the sparse matrix
kx = 3*mnx;
ky = 3*mny;
kz = 3*mnz;
Mx = sparse(iix(1:kx), jjx(1:kx), sx(1:kx), (Nr+2)*(Ntheta+2)*(Nphi+2), (Nr+2)*(Ntheta+2)*(Nphi+2));
My = sparse(iiy(1:ky), jjy(1:ky), sy(1:ky), (Nr+2)*(Ntheta+2)*(Nphi+2), (Nr+2)*(Ntheta+2)*(Nphi+2));
Mz = sparse(iiz(1:kz), jjz(1:kz), sz(1:kz), (Nr+2)*(Ntheta+2)*(Nphi+2), (Nr+2)*(Ntheta+2)*(Nphi+2));
M = Mx + My + Mz;
end
