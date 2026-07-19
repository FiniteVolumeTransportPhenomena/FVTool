function [RHS, RHSx, RHSy, RHSz] = ...
    convectionTvdRHSSpherical3D(u, phi, FL)
% This function uses the TVD scheme to discretize the explicit correction
% (RHS) of a 3D convection term \grad . (u \phi) on a spherical
% (r, theta, phi) mesh, where u is a face vector and FL a flux limiter.
%
% SYNOPSIS:
%   [RHS, RHSx, RHSy, RHSz] = convectionTvdRHSSpherical3D(u, phi, FL)
%
% SEE ALSO: convectionTvdRHSCylindrical3D

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
sinf_n = sin(thetaf(:,2:Ntheta+1,:));
sinf_s = sin(thetaf(:,1:Ntheta,:));
dx=repmat(0.5*(u.domain.cellsize.x(1:end-1)+u.domain.cellsize.x(2:end)), 1, Ntheta, Nphi);
dy=repmat(0.5*(u.domain.cellsize.y(1:end-1)+u.domain.cellsize.y(2:end))', Nr, 1, Nphi);
dz=zeros(1, 1, Nphi+1);
dz(1,1,:)=0.5*(u.domain.cellsize.z(1:end-1)+u.domain.cellsize.z(2:end));
dz=repmat(dz, Nr, Ntheta, 1);
psiX_p = zeros(Nr+1,Ntheta,Nphi);
psiX_m = zeros(Nr+1,Ntheta,Nphi);
psiY_p = zeros(Nr,Ntheta+1,Nphi);
psiY_m = zeros(Nr,Ntheta+1,Nphi);
psiZ_p = zeros(Nr,Ntheta,Nphi+1);
psiZ_m = zeros(Nr,Ntheta,Nphi+1);

% define the vectors to stores the sparse matrix data
mnx = Nr*Ntheta*Nphi;	mny = Nr*Ntheta*Nphi;   mnz = Nr*Ntheta*Nphi;

% extract the velocity data
ux = u.xvalue;
uy = u.yvalue;
uz = u.zvalue;

% calculate the upstream to downstream gradient ratios for u>0 (+ ratio)
% x direction
dphiX_p = (phi.value(2:Nr+2, 2:Ntheta+1, 2:Nphi+1)-phi.value(1:Nr+1, 2:Ntheta+1, 2:Nphi+1))./dx;
rX_p = dphiX_p(1:end-1,:,:)./fsign(dphiX_p(2:end,:,:));
psiX_p(2:Nr+1,:,:) = 0.5*FL(rX_p).* ...
    (phi.value(3:Nr+2,2:Ntheta+1,2:Nphi+1)-phi.value(2:Nr+1,2:Ntheta+1,2:Nphi+1));
psiX_p(1,:,:) = 0; % left boundary
% y direction
dphiY_p = (phi.value(2:Nr+1, 2:Ntheta+2, 2:Nphi+1)-phi.value(2:Nr+1, 1:Ntheta+1, 2:Nphi+1))./dy;
rY_p = dphiY_p(:,1:end-1,:)./fsign(dphiY_p(:,2:end,:));
psiY_p(:,2:Ntheta+1,:) = 0.5*FL(rY_p).* ...
    (phi.value(2:Nr+1,3:Ntheta+2,2:Nphi+1)-phi.value(2:Nr+1, 2:Ntheta+1,2:Nphi+1));
psiY_p(:,1,:) = 0; % Bottom boundary
% z direction
dphiZ_p = (phi.value(2:Nr+1, 2:Ntheta+1, 2:Nphi+2)-phi.value(2:Nr+1, 2:Ntheta+1, 1:Nphi+1))./dz;
rZ_p = dphiZ_p(:,:,1:end-1)./fsign(dphiZ_p(:,:,2:end));
psiZ_p(:,:,2:Nphi+1) = 0.5*FL(rZ_p).* ...
    (phi.value(2:Nr+1,2:Ntheta+1,3:Nphi+2)-phi.value(2:Nr+1,2:Ntheta+1,2:Nphi+1));
psiZ_p(:,:,1) = 0; % Back boundary

% calculate the upstream to downstream gradient ratios for u<0 (- ratio)
% x direction
rX_m = dphiX_p(2:end,:,:)./fsign(dphiX_p(1:end-1,:,:));
psiX_m(1:Nr,:,:) = 0.5*FL(rX_m).* ...
    (phi.value(1:Nr, 2:Ntheta+1, 2:Nphi+1)-phi.value(2:Nr+1, 2:Ntheta+1, 2:Nphi+1));
psiX_m(Nr+1,:,:) = 0; % right boundary
% y direction
rY_m = dphiY_p(:,2:end,:)./fsign(dphiY_p(:,1:end-1,:));
psiY_m(:,1:Ntheta,:) = 0.5*FL(rY_m).* ...
    (phi.value(2:Nr+1,1:Ntheta,2:Nphi+1)-phi.value(2:Nr+1,2:Ntheta+1,2:Nphi+1));
psiY_m(:,Ntheta+1,:) = 0; % top boundary
% z direction
rZ_m = dphiZ_p(:,:,2:end)./fsign(dphiZ_p(:,:,1:end-1));
psiZ_m(:,:,1:Nphi) = 0.5*FL(rZ_m).* ...
    (phi.value(2:Nr+1,2:Ntheta+1,1:Nphi)-phi.value(2:Nr+1,2:Ntheta+1,2:Nphi+1));
psiZ_m(:,:,Nphi+1) = 0; % front boundary
% reassign the east, west, north, south, front, back velocity vectors
ue = ux(2:Nr+1,:,:);		uw = ux(1:Nr,:,:);
vn = uy(:,2:Ntheta+1,:);     vs = uy(:,1:Ntheta,:);
wf = uz(:,:,2:Nphi+1);     wb = uz(:,:,1:Nphi);
re = rf(2:Nr+1,:,:);         rw = rf(1:Nr,:,:);

% find the velocity direction for the upwind scheme
ue_min = min(ue,0);	ue_max = max(ue,0);
uw_min = min(uw,0);	uw_max = max(uw,0);
vn_min = min(vn,0);	vn_max = max(vn,0);
vs_min = min(vs,0);	vs_max = max(vs,0);
wf_min = min(wf,0);	wf_max = max(wf,0);
wb_min = min(wb,0);	wb_max = max(wb,0);

% build the sparse matrix based on the numbering system
rowx_index = reshape(G(2:Nr+1,2:Ntheta+1,2:Nphi+1),mnx,1); % main diagonal x
rowy_index = reshape(G(2:Nr+1,2:Ntheta+1,2:Nphi+1),mny,1); % main diagonal y
rowz_index = reshape(G(2:Nr+1,2:Ntheta+1,2:Nphi+1),mnz,1); % main diagonal z

% calculate the TVD correction term with the spherical metric factors
div_x = -(1./(DRp.*rp.^2)).*(re.^2.*(ue_max.*psiX_p(2:Nr+1,:,:)+ue_min.*psiX_m(2:Nr+1,:,:))- ...
              rw.^2.*(uw_max.*psiX_p(1:Nr,:,:)+uw_min.*psiX_m(1:Nr,:,:)));
div_y = -(1./(DTHETAp.*rp.*sinp)).*(sinf_n.*(vn_max.*psiY_p(:,2:Ntheta+1,:)+vn_min.*psiY_m(:,2:Ntheta+1,:))- ...
              sinf_s.*(vs_max.*psiY_p(:,1:Ntheta,:)+vs_min.*psiY_m(:,1:Ntheta,:)));
div_z = -(1./(DPHIp.*rp.*sinp)).*((wf_max.*psiZ_p(:,:,2:Nphi+1)+wf_min.*psiZ_m(:,:,2:Nphi+1))- ...
              (wb_max.*psiZ_p(:,:,1:Nphi)+wb_min.*psiZ_m(:,:,1:Nphi)));

% define the RHS Vector
RHS = zeros((Nr+2)*(Ntheta+2)*(Nphi+2),1);
RHSx = zeros((Nr+2)*(Ntheta+2)*(Nphi+2),1);
RHSy = zeros((Nr+2)*(Ntheta+2)*(Nphi+2),1);
RHSz = zeros((Nr+2)*(Ntheta+2)*(Nphi+2),1);

% assign the values of the RHS vector
row_index = rowx_index;
RHS(row_index) = reshape(div_x+div_y+div_z,Nr*Ntheta*Nphi,1);
RHSx(rowx_index) = reshape(div_x,Nr*Ntheta*Nphi,1);
RHSy(rowy_index) = reshape(div_y,Nr*Ntheta*Nphi,1);
RHSz(rowz_index) = reshape(div_z,Nr*Ntheta*Nphi,1);

end

function phi_out = fsign(phi_in)
% This function checks the value of phi_in and assigns an eps value to the
% elements that are less than or equal to zero, while keeping the signs of
% the nonzero elements
    phi_out = (abs(phi_in)>=eps).*phi_in+eps*(phi_in==0)+eps*(abs(phi_in)<eps).*sign(phi_in);
end
