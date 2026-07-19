function [M, Mx, My, Mz] = convectionUpwindTermSpherical3D(u, varargin)
% This function uses the upwind scheme to discretize a 3D convection term
% in the form \grad . (u \phi) on a spherical (r, theta, phi) mesh, where
% u is a face vector. An optional second face vector selects the upwind
% direction (defaults to u itself).
%
% SYNOPSIS:
%   [M, Mx, My, Mz] = convectionUpwindTermSpherical3D(u)
%   [M, Mx, My, Mz] = convectionUpwindTermSpherical3D(u, u_upwind)
%
% SEE ALSO: convectionUpwindTermCylindrical3D

% extract data from the mesh structure
Nr = u.domain.dims(1);
Ntheta = u.domain.dims(2);
Nphi = u.domain.dims(3);
G=reshape((1:(Nr+2)*(Ntheta+2)*(Nphi+2)), Nr+2, Ntheta+2, Nphi+2);
DRp = repmat(u.domain.cellsize.x(2:end-1), 1, Ntheta, Nphi);
DTHETAp = repmat(u.domain.cellsize.y(2:end-1)', Nr, 1, Nphi);
DPHI = zeros(1,1,Nphi+2);
DPHI(1,1,:) = u.domain.cellsize.z;
DPHIp=repmat(DPHI(1,1,2:end-1), Nr, Ntheta, 1);
rp = repmat(u.domain.cellcenters.x, 1, Ntheta, Nphi);
rf = repmat(u.domain.facecenters.x, 1, Ntheta, Nphi);
thetap = repmat(u.domain.cellcenters.y', Nr, 1, Nphi);
thetaf = repmat(u.domain.facecenters.y', Nr, 1, Nphi);
sinp = sin(thetap);

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

% extract the velocity data; use u_upwind for the flow direction
ux = u.xvalue;   uy = u.yvalue;   uz = u.zvalue;
if nargin>1
    ux_up = varargin{1}.xvalue; uy_up = varargin{1}.yvalue; uz_up = varargin{1}.zvalue;
else
    ux_up = ux; uy_up = uy; uz_up = uz;
end

% reassign the east, west, north, south, front, back velocity vectors
re = rf(2:Nr+1,:,:);         rw = rf(1:Nr,:,:);
sinf_n = sin(thetaf(:,2:Ntheta+1,:));   sinf_s = sin(thetaf(:,1:Ntheta,:));

% find the velocity direction for the upwind scheme (from u_upwind)
ue_min = min(ux_up(2:Nr+1,:,:),0);	ue_max = max(ux_up(2:Nr+1,:,:),0);
uw_min = min(ux_up(1:Nr,:,:),0);	uw_max = max(ux_up(1:Nr,:,:),0);
vn_min = min(uy_up(:,2:Ntheta+1,:),0);	vn_max = max(uy_up(:,2:Ntheta+1,:),0);
vs_min = min(uy_up(:,1:Ntheta,:),0);	vs_max = max(uy_up(:,1:Ntheta,:),0);
wf_min = min(uz_up(:,:,2:Nphi+1),0);	wf_max = max(uz_up(:,:,2:Nphi+1),0);
wb_min = min(uz_up(:,:,1:Nphi),0);	wb_max = max(uz_up(:,:,1:Nphi),0);

% calculate the coefficients for the internal cells
AE = re.^2.*ue_min./(DRp.*rp.^2);
AW = -rw.^2.*uw_max./(DRp.*rp.^2);
AN = vn_min.*sinf_n./(DTHETAp.*rp.*sinp);
AS = -vs_max.*sinf_s./(DTHETAp.*rp.*sinp);
AF = wf_min./(DPHIp.*rp.*sinp);
AB = -wb_max./(DPHIp.*rp.*sinp);
APx = (re.^2.*ue_max-rw.^2.*uw_min)./(DRp.*rp.^2);
APy = (sinf_n.*vn_max-sinf_s.*vs_min)./(DTHETAp.*rp.*sinp);
APz = (wf_max-wb_min)./(DPHIp.*rp.*sinp);

% Also correct for the boundary cells (not the ghost cells)
% Left boundary:
APx(1,:,:) = APx(1,:,:)-rw(1,:,:).^2.*uw_max(1,:,:)./(2*DRp(1,:,:).*rp(1,:,:).^2);   AW(1,:,:) = AW(1,:,:)/2;
% Right boundary:
AE(end,:,:) = AE(end,:,:)/2;    APx(end,:,:) = APx(end,:,:)+re(end,:,:).^2.*ue_min(end,:,:)./(2*DRp(end,:,:).*rp(end,:,:).^2);
% Bottom boundary:
APy(:,1,:) = APy(:,1,:)-vs_max(:,1,:).*sinf_s(:,1,:)./(2*DTHETAp(:,1,:).*rp(:,1,:).*sinp(:,1,:));   AS(:,1,:) = AS(:,1,:)/2;
% Top boundary:
AN(:,end,:) = AN(:,end,:)/2;    APy(:,end,:) = APy(:,end,:)+vn_min(:,end,:).*sinf_n(:,end,:)./(2*DTHETAp(:,end,:).*rp(:,end,:).*sinp(:,end,:));
% Back boundary:
APz(:,:,1) = APz(:,:,1)-wb_max(:,:,1)./(2*DPHIp(:,:,1).*rp(:,:,1).*sinp(:,:,1));   AB(:,:,1) = AB(:,:,1)/2;
% Front boundary:
AF(:,:,end) = AF(:,:,end)/2;    APz(:,:,end) = APz(:,:,end)+wf_min(:,:,end)./(2*DPHIp(:,:,end).*rp(:,:,end).*sinp(:,:,end));

AE = reshape(AE,mnx,1);
AW = reshape(AW,mnx,1);
AN = reshape(AN,mny,1);
AS = reshape(AS,mny,1);
AF = reshape(AF,mnz,1);
AB = reshape(AB,mnz,1);
APx = reshape(APx,mnx,1);
APy = reshape(APy,mny,1);
APz = reshape(APz,mnz,1);

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
