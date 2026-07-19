function cellvol = cellVolume(meshvar)
% cellvol = cellVolume(meshvar)
% returns the volume of each cell as a cell variable
% SYNOPSIS:
%   cellvol = cellVolume(meshvar)
%
% PARAMETERS:
%   MeshStructure: a mesh structure created by buildMesh* functions
%
% RETURNS:
%   cellvol: a (1D, 2D, or 3D) matrix depending on the mesh size
%
% EXAMPLE:
%   m = createMesh2D(3,4, 1.0, 2.0); % creates a mesh
%   cell_vol=cellVolume(m);
%
% SEE ALSO:
%     createFaceVariable, createBC, buildMesh1D,
%     buildMesh2D, buildMesh3D,
%     buildMeshCylindrical1D, buildMeshCylindrical2D,
%     cellBoundary, combineBC
%
% Written by Ali A. Eftekhari

% check the size of the variable and the mesh dimension
dim = meshvar.dimension;
BC = createBC(meshvar);
switch geometryTag(meshvar)
    case '1D'
        c=meshvar.cellsize.x(2:end-1);
    case 'Cylindrical1D'
        c=2.0*pi()*meshvar.cellsize.x(2:end-1).*meshvar.cellcenters.x;
    case 'Spherical1D'
        c=4.0*pi()*meshvar.cellcenters.x.^2.*meshvar.cellsize.x(2:end-1);
    case '2D'
        c=meshvar.cellsize.x(2:end-1)*meshvar.cellsize.y(2:end-1)';
    case 'Cylindrical2D' % cylindrical
        c=2.0*pi()*meshvar.cellcenters.x.*meshvar.cellsize.x(2:end-1)*meshvar.cellsize.y(2:end-1)';
    case 'Radial2D' % radial
        c=meshvar.cellcenters.x.*meshvar.cellsize.x(2:end-1)*meshvar.cellsize.y(2:end-1)';
    case '3D'
        Nx = meshvar.dims(1);
        Ny = meshvar.dims(2);
        Nz = meshvar.dims(3);
        DXp = repmat(meshvar.cellsize.x(2:end-1), 1, Ny, Nz);
        DYp = repmat(meshvar.cellsize.y(2:end-1)', Nx, 1, Nz);
        DZ = zeros(1,1,Nz+2);
        DZ(1,1,:) = meshvar.cellsize.z;
        DZp=repmat(DZ(1,1,2:end-1), Nx, Ny, 1);
        c=DXp.*DYp.*DZp;
    case 'Cylindrical3D'
        N = meshvar.dims;
        Nr = N(1); Ntetta=N(2); Nz = N(3);
        rp = repmat(meshvar.cellcenters.x, 1, Ntetta, Nz);
        DRp = repmat(meshvar.cellsize.x(2:end-1), 1, Ntetta, Nz);
        DTHETAp = repmat(meshvar.cellsize.y(2:end-1)', Nr, 1, Nz);
        DZ = zeros(1,1,Nz+2);
        DZ(1,1,:) = meshvar.cellsize.z;
        DZp=repmat(DZ(1,1,2:end-1), Nr, Ntetta, 1);
        c=rp.*DRp.*DTHETAp.*DZp;
    case 'Spherical3D'
        N = meshvar.dims;
        Nr = N(1); Ntheta=N(2); Nphi = N(3);
        rp = repmat(meshvar.cellcenters.x, 1, Ntheta, Nphi);
        THETAp = repmat(meshvar.cellcenters.y', Nr, 1, Nphi);
        DRp = repmat(meshvar.cellsize.x(2:end-1), 1, Ntheta, Nphi);
        DTHETAp = repmat(meshvar.cellsize.y(2:end-1)', Nr, 1, Nphi);
        DPHI = zeros(1,1,Nphi+2);
        DPHI(1,1,:) = meshvar.cellsize.z;
        DPHIp = repmat(DPHI(1,1,2:end-1), Nr, Ntheta, 1);
        c=rp.^2.*sin(THETAp).*DRp.*DTHETAp.*DPHIp;
    otherwise
        error('FVTool:unsupportedGeometry', ...
            'cellVolume: no implementation for %s', geometryTag(meshvar));
end
cellvol= createCellVariable(meshvar, c, BC);
