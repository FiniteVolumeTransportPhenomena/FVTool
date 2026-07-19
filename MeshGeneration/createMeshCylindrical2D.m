function MS = createMeshCylindrical2D(varargin)
% MeshStructure = createMeshCylindrical2D(Nr, Ny, Lr, Ly)
% MeshStructure = createMeshCylindrical2D(facelocationR, facelocationY)
% creates a uniform 2D mesh on a cylindrical coordinate:
% Nr is the number of cells in r (radial) direction
% Ny is the number of cells in y (vertical) direction
% Lr is the domain length in r direction
% Ly is the domain length in y direction
%
% SYNOPSIS:
%   MeshStructure = createMeshCylindrical2D(Nr, Ny, Lr, Ly)
%
% PARAMETERS:
%   Nr: number of cells in the x direction
%   Lr: domain length in x direction
%   Ny: number of cells in the y direction
%   Ly: domain length in y direction
%
% RETURNS:
%   MeshStructure.
%                 dimensions=2.5 (2D problem, cylindrical coordinate)
%                 numbering: shows the indexes of cellsn from left to right
%                 and top to bottom
%                 cellsize: r and y elements of the cell size =[Lr/Nr,
%                 Ly/Ny]
%                 cellcenters.x: location of each cell in the r direction
%                 cellcenters.y: location of each cell in the y direction
%                 facecenters.x: location of interface between cells in the
%                 r direction
%                 facecenters.y: location of interface between cells in the
%                 y direction
%                 numberofcells: [Nr, Ny]
%
%
% EXAMPLE:
%   Nr = 5;
%   Ny = 7;
%   R = 10;
%   Ly = 20;
%   m = createMeshCylindrical2D(Nx, Ny, Lx, Ly);
%   [X, Y] = ndgrid(m.cellcenters.x, m.cellcenters.y);
%   [Xf,Yf]=ndgrid(m.facecenters.x, m.facecenters.y);
%   plot(X, Y, 'or', ...
%        Xf, Yf, '-b', Xf', Yf', '-b');
%
% SEE ALSO:
%     createMesh1D, createMesh3D, createMeshCylindrical1D, ...
%     createMeshCylindrical2D, createCellVariable, createFaceVariable

% Written by Ali A. Eftekhari
% See the license file

if nargin==4
  % uniform 1D mesh
  Nx=varargin{1};
  Ny=varargin{2};
  Width=varargin{3};
  Height=varargin{4};
  % cell size is dx
  dx = Width/Nx;
  dy = Height/Ny;
  G=reshape(1:(Nx+2)*(Ny+2), Nx+2, Ny+2);
  CellSize.x= dx*ones(Nx+2,1);
  CellSize.y= dy*ones(Ny+2,1);
  CellSize.z= [0.0];
  CellLocation.x= [1:Nx]'*dx-dx/2;
  CellLocation.y= [1:Ny]'*dy-dy/2;
  CellLocation.z= [0.0];
  FaceLocation.x= [0:Nx]'*dx;
  FaceLocation.y= [0:Ny]'*dy;
  FaceLocation.z= [0.0];
elseif nargin==2
  % nonuniform 1D mesh
  facelocationX=varargin{1};
  facelocationY=varargin{2};
  facelocationX=facelocationX(:);
  facelocationY=facelocationY(:);
  Nx = length(facelocationX)-1;
  Ny = length(facelocationY)-1;
  G=reshape(1:(Nx+2)*(Ny+2), Nx+2, Ny+2);
  CellSize.x= [facelocationX(2)-facelocationX(1); ...
    facelocationX(2:end)-facelocationX(1:end-1); ...
    facelocationX(end)-facelocationX(end-1)];
  CellSize.y= [facelocationY(2)-facelocationY(1); ...
    facelocationY(2:end)-facelocationY(1:end-1); ...
    facelocationY(end)-facelocationY(end-1)];
  CellSize.z= [0.0];
  CellLocation.x= 0.5*(facelocationX(2:end)+facelocationX(1:end-1));
  CellLocation.y= 0.5*(facelocationY(2:end)+facelocationY(1:end-1));
  CellLocation.z= [0.0];
  FaceLocation.x= facelocationX;
  FaceLocation.y= facelocationY;
  FaceLocation.z= [0.0];
end
c=G([1,end], [1,end]);
MS=MeshCylindrical2D( ...
  [Nx,Ny], ...
  CellSize, ...
  CellLocation, ...
  FaceLocation, ...
  c(:), ...
  [1]);
