function MS = createMesh1D(varargin)
% MeshStructure = createMesh1D(Nx, Width)
% MeshStructure = createMesh1D(facelocationX)
% creates a uniform 1D mesh:
% Nx is the number of cells in x (horizontal) direction
% Width is the domain length in x direction
%
% SYNOPSIS:
%   MeshStructure = createMesh1D(Nx, Width)
%
% PARAMETERS:
%   Nx: number of cells in the domain
%   Width: domain length
%
% RETURNS:
%   MeshStructure.
%                 dimensions=1 (1D problem)
%                 numbering: shows the indexes of cellsn from left to right
%                 dx: cell size (=Width/Nx)
%                 cellcenters.x: location of each cell in the x direction
%                 facecenters.x: location of interface between cells in the
%                 x direction
%                 numberofcells: Nx
%
%
% EXAMPLE:
%   L = 1.0; % length of the domain
%   Nx = 10; % number of cells in the domain
%   m = createMesh1D(Nx, L);
%   plot(m.cellcenters.x, ones(size(m.cellcenters.x)), 'o', ...
%        m.facecenters.x, ones(size(m.facecenters.x)), '-+');
%   legend('cell centers', 'face centers');
%
% SEE ALSO:
%     createMesh2D, createMesh3D, createMeshCylindrical1D, ...
%     createMeshCylindrical2D, createCellVariable, createFaceVariable

% Written by Ali A. Eftekhari
% See the license file

if nargin==2
  % uniform 1D mesh
  Nx=varargin{1};
  Width=varargin{2};
  % cell size is dx
  dx = Width/Nx;
  CellSize.x= dx*ones(Nx+2,1);
  CellSize.y= [0.0];
  CellSize.z= [0.0];
  CellLocation.x= [1:Nx]'*dx-dx/2;
  CellLocation.y= [0.0];
  CellLocation.z= [0.0];
  FaceLocation.x= [0:Nx]'*dx;
  FaceLocation.y= [0.0];
  FaceLocation.z= [0.0];
elseif nargin==1
  % nonuniform 1D mesh
  facelocationX=varargin{1};
  n=size(facelocationX);
  if n(1)==1
      facelocationX=facelocationX';
  end
  Nx = length(facelocationX)-1;
  CellSize.x= [facelocationX(2)-facelocationX(1); ...
    facelocationX(2:end)-facelocationX(1:end-1); ...
    facelocationX(end)-facelocationX(end-1)];
  CellSize.y= [0.0];
  CellSize.z= [0.0];
  CellLocation.x= 0.5*(facelocationX(2:end)+facelocationX(1:end-1));
  CellLocation.y= [0.0];
  CellLocation.z= [0.0];
  FaceLocation.x= facelocationX;
  FaceLocation.y= [0.0];
  FaceLocation.z= [0.0];
end

MS=Mesh1D( ...
  [Nx,1], ...
  CellSize, ...
  CellLocation, ...
  FaceLocation, ...
  [1], ...
  [1]);
