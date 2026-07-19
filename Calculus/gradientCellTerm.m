function cellGrad = gradientCellTerm(phi)
% this function calculates the gradient of a variable in x direction, in
% the cell center. It needs a cell variable as an input.
%
% SYNOPSIS:
%   cellGrad = gradientCellTerm(phi)
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


phi_face=linearMean(phi);
tag = geometryTag(phi.domain);
if any(strcmp(tag, {'1D', 'Cylindrical1D', 'Spherical1D'}))
	Nx = phi.domain.dims(1);
	DX = phi.domain.cellsize.x(2:end-1);
	xvalue = (phi_face.xvalue(2:Nx+1)-phi_face.xvalue(1:Nx))./DX;
	yvalue=[];
	zvalue=[];
	cellGrad=FaceVariable(phi.domain, xvalue, yvalue, zvalue);
elseif any(strcmp(tag, {'2D', 'Cylindrical2D'}))
	Nx = phi.domain.dims(1);
	Ny = phi.domain.dims(2);
	DX = repmat(phi.domain.cellsize.x(2:end-1), 1, Ny);
	DY = repmat(phi.domain.cellsize.y(2:end-1)', Nx, 1);
	xvalue = (phi_face.xvalue(2:Nx+1,:)-phi_face.xvalue(1:Nx,:))./DX;
	yvalue = (phi_face.yvalue(:,2:Ny+1)-phi_face.yvalue(:,1:Ny))./DY;
	zvalue=[];
	cellGrad=FaceVariable(phi.domain, xvalue, yvalue, zvalue);
elseif strcmp(tag, 'Radial2D')
	Nr = phi.domain.dims(1);
	Ntheta = phi.domain.dims(2);
	DR = repmat(phi.domain.cellsize.x(2:end-1), 1, Ntheta);
	DTHETA = repmat(phi.domain.cellsize.y(2:end-1)', Nr, 1);
	rp = repmat(phi.domain.cellcenters.x, 1, Ntheta);
	xvalue = (phi_face.xvalue(2:Nr+1,:)-phi_face.xvalue(1:Nr,:))./DR;
	yvalue = (phi_face.yvalue(:,2:Ntheta+1)-phi_face.yvalue(:,1:Ntheta))./(DTHETA.*rp);
	zvalue=[];
	cellGrad=FaceVariable(phi.domain, xvalue, yvalue, zvalue);
elseif strcmp(tag, '3D')
	Nx = phi.domain.dims(1);
	Ny = phi.domain.dims(2);
	Nz = phi.domain.dims(3);
	DX = repmat(phi.domain.cellsize.x(2:end-1), 1, Ny, Nz);
	DY = repmat(phi.domain.cellsize.y(2:end-1)', Nx, 1, Nz);
	DZ = zeros(1,1,Nz);
	DZ(1,1,:) = phi.domain.cellsize.z(2:end-1);
	DZ=repmat(DZ, Nx, Ny, 1);
	xvalue = (phi_face.xvalue(2:Nx+1,:,:)-phi_face.xvalue(1:Nx,:,:))./DX;
	yvalue = (phi_face.yvalue(:,2:Ny+1,:)-phi_face.yvalue(:,1:Ny,:))./DY;
	zvalue = (phi_face.zvalue(:,:,2:Nz+1)-phi_face.zvalue(:,:,1:Nz))./DZ;
	cellGrad=FaceVariable(phi.domain, xvalue, yvalue, zvalue);
elseif strcmp(tag, 'Cylindrical3D')
	Nr = phi.domain.dims(1);
	Ntheta = phi.domain.dims(2);
	Nz = phi.domain.dims(3);
	DR = repmat(phi.domain.cellsize.x(2:end-1), 1, Ntheta, Nz);
	DTHETA = repmat(phi.domain.cellsize.y(2:end-1)', Nr, 1, Nz);
	DZ = zeros(1,1,Nz);
	DZ(1,1,:) = phi.domain.cellsize.z(2:end-1);
	DZ = repmat(DZ, Nr, Ntheta, 1);
	rp = repmat(phi.domain.cellcenters.x, 1, Ntheta, Nz);
	xvalue = (phi_face.xvalue(2:Nr+1,:,:)-phi_face.xvalue(1:Nr,:,:))./DR;
	yvalue = (phi_face.yvalue(:,2:Ntheta+1,:)-phi_face.yvalue(:,1:Ntheta,:))./(DTHETA.*rp);
	zvalue = (phi_face.zvalue(:,:,2:Nz+1)-phi_face.zvalue(:,:,1:Nz))./DZ;
	cellGrad=FaceVariable(phi.domain, xvalue, yvalue, zvalue);
elseif strcmp(tag, 'Spherical3D')
	Nr = phi.domain.dims(1);
	Ntheta = phi.domain.dims(2);
	Nphi = phi.domain.dims(3);
	DR = repmat(phi.domain.cellsize.x(2:end-1), 1, Ntheta, Nphi);
	DTHETA = repmat(phi.domain.cellsize.y(2:end-1)', Nr, 1, Nphi);
	DPHI = zeros(1,1,Nphi);
	DPHI(1,1,:) = phi.domain.cellsize.z(2:end-1);
	DPHI = repmat(DPHI, Nr, Ntheta, 1);
	rp = repmat(phi.domain.cellcenters.x, 1, Ntheta, Nphi);
	thetap = repmat(phi.domain.cellcenters.y', Nr, 1, Nphi);
	xvalue = (phi_face.xvalue(2:Nr+1,:,:)-phi_face.xvalue(1:Nr,:,:))./DR;
	yvalue = (phi_face.yvalue(:,2:Ntheta+1,:)-phi_face.yvalue(:,1:Ntheta,:))./(DTHETA.*rp);
	zvalue = (phi_face.zvalue(:,:,2:Nphi+1)-phi_face.zvalue(:,:,1:Nphi))./(DPHI.*rp.*sin(thetap));
	cellGrad=FaceVariable(phi.domain, xvalue, yvalue, zvalue);
else
	error('FVTool:unsupportedGeometry', ...
		'gradientCellTerm: no implementation for %s', tag);
end
