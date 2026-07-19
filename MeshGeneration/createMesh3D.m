function MS = createMesh3D(varargin)
% MeshStructure = createMesh3D(Nx, Ny, Nz, Width, Height, Depth)
% creates a uniform 3D mesh:
% Nx is the number of cells in x (horizontal) direction
% Ny is the number of cells in y (vertical) direction
% Nz is the number of cells in z (perpendicular) direction
% Lx is the domain length in x direction
% Ly is the domain length in y direction
% Lz is the domain length in z direction
%
% SYNOPSIS:
%   MeshStructure = createMesh3D(Nx, Ny, Nz, Lx, Ly, Lz)
%
% PARAMETERS:
%   Nx: number of cells in the x direction
%   Lx: domain length in x direction
%   Ny: number of cells in the y direction
%   Ly: domain length in y direction
%   Nz: number of cells in the z direction
%   Lz: domain length in z direction
%
% RETURNS:
%   MeshStructure.
%                 dimensions=3 (3D problem)
%                 numbering: shows the indexes of cellsn from left to right
%                 and top to bottom and back to front
%                 cellsize: x, y, and z elements of the cell size =[Lx/Nx,
%                 Ly/Ny, Lz/Nz]
%                 cellcenters.x: location of each cell in the x direction
%                 cellcenters.y: location of each cell in the y direction
%                 cellcenters.z: location of each cell in the z direction
%                 facecenters.x: location of interface between cells in the
%                 x direction
%                 facecenters.y: location of interface between cells in the
%                 y direction
%                 facecenters.z: location of interface between cells in the
%                 z direction
%                 numberofcells: [Nx, Ny, Nz]
%
%
% EXAMPLE:
%   Nx = 2;
%   Lx = 1.0;
%   Ny = 3;
%   Ly = 2.0;
%   Nz = 4;
%   Lz = 3.0;
%   m = createMesh3D(Nx, Ny, Nz, Lx, Ly, Lz);
%   [X, Y, Z] = ndgrid(m.cellcenters.x, m.cellcenters.y, m.cellcenters.z);
%   [Xf, Yf, Zf] = ndgrid(m.facecenters.x, m.facecenters.y, m.facecenters.z);
%   plot3(X(:), Y(:), Z(:), 'or')
%   hold on;
%   plot3(Xf(:), Yf(:), Zf(:), '+b')
%   legend('cell centers', 'cell corners');
%
% SEE ALSO:
%     createMesh1D, createMesh2D, createMeshCylindrical1D, ...
%     createMeshCylindrical2D, createCellVariable, createFaceVariable

% Written by Ali A. Eftekhari
% See the license file

if nargin==6
  % uniform 1D mesh
  Nx=varargin{1};
  Ny=varargin{2};
  Nz=varargin{3};
  Width=varargin{4};
  Height=varargin{5};
  Depth=varargin{6};
  % cell size is dx
  dx = Width/Nx;
  dy = Height/Ny;
  dz = Depth/Nz;
  G=reshape(1:(Nx+2)*(Ny+2)*(Nz+2), Nx+2, Ny+2, Nz+2);
  CellSize.x= dx*ones(Nx+2,1);
  CellSize.y= dy*ones(Ny+2,1);
  CellSize.z= dz*ones(Nz+2,1);
  CellLocation.x= [1:Nx]'*dx-dx/2;
  CellLocation.y= [1:Ny]'*dy-dy/2;
  CellLocation.z= [1:Nz]'*dz-dz/2;
  FaceLocation.x= [0:Nx]'*dx;
  FaceLocation.y= [0:Ny]'*dy;
  FaceLocation.z= [0:Nz]'*dz;
elseif nargin==3
  % nonuniform 1D mesh
  facelocationX=varargin{1};
  facelocationY=varargin{2};
  facelocationZ=varargin{3};
  facelocationX=facelocationX(:);
  facelocationY=facelocationY(:);
  facelocationZ=facelocationZ(:);
  Nx = length(facelocationX)-1;
  Ny = length(facelocationY)-1;
  Nz = length(facelocationZ)-1;
  G=reshape(1:(Nx+2)*(Ny+2)*(Nz+2), Nx+2, Ny+2, Nz+2);
  CellSize.x= [facelocationX(2)-facelocationX(1); ...
    facelocationX(2:end)-facelocationX(1:end-1); ...
    facelocationX(end)-facelocationX(end-1)];
  CellSize.y= [facelocationY(2)-facelocationY(1); ...
    facelocationY(2:end)-facelocationY(1:end-1); ...
    facelocationY(end)-facelocationY(end-1)];
  CellSize.z= [facelocationZ(2)-facelocationZ(1); ...
    facelocationZ(2:end)-facelocationZ(1:end-1); ...
    facelocationZ(end)-facelocationZ(end-1)];
  CellLocation.x= 0.5*(facelocationX(2:end)+facelocationX(1:end-1));
  CellLocation.y= 0.5*(facelocationY(2:end)+facelocationY(1:end-1));
  CellLocation.z= 0.5*(facelocationZ(2:end)+facelocationZ(1:end-1));
  FaceLocation.x= facelocationX;
  FaceLocation.y= facelocationY;
  FaceLocation.z= facelocationZ;
end
c=G([1,end], [1,end], [1, end]);
e1=G([1, end], [1, end], 2:Nz+1);
e2=G([1, end], 2:Ny+1, [1, end]);
e3=G(2:Nx+1, [1, end], [1, end]);
MS=Mesh3D( ...
  [Nx,Ny, Nz], ...
  CellSize, ...
  CellLocation, ...
  FaceLocation, ...
  c(:), ...
  [e1(:); e2(:); e3(:)]);
