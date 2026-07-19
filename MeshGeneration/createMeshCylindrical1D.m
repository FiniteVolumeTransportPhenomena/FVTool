function MS = createMeshCylindrical1D(varargin)
% MeshStructure = createMeshCylindrical1D(Nr, Lr)
% MeshStructure = createMeshCylindrical1D(facelocationR)
% creates a uniform 1D mesh (1D axial symmetry):
% Nx is the number of cells in r (radial) direction
% Width is the domain radius in r direction
%
% SYNOPSIS:
%   MeshStructure = createMeshCylindrical1D(Nr, Lr)
%
% PARAMETERS:
%   Nr: number of cells in the domain
%   Lr: domain radius
%
% RETURNS:
%   MeshStructure.
%                 dimensions=1.5 (1D axial symmetry)
%                 numbering: shows the indexes of cellsn from left to right
%                 dx: cell size (=Lr/Nx)
%                 cellcenters.x: location of each cell in the x direction
%                 facecenters.x: location of interface between cells in the
%                 r direction
%                 numberofcells: Nr
%
%
% EXAMPLE:
%   R = 1.0; % length of the domain
%   Nr = 10; % number of cells in the domain
%   m = createMeshCylindrical1D(Nr, R);
%   plot(m.cellcenters.x, ones(size(m.cellcenters.x)), 'o', ...
%        m.facecenters.x, ones(size(m.facecenters.x)), '-+');
%   legend('cell centers', 'face centers');
%
% SEE ALSO:
%     createMesh1D, createMesh2D, createMesh3D, ...
%     createMeshCylindrical2D, createCellVariable, createFaceVariable

% Written by Ali A. Eftekhari
% See the license file

if nargin==2
  % uniform 1D mesh
  Nr=varargin{1};
  R=varargin{2};
  % cell size is dx
  dx = R/Nr;
  CellSize.x= dx*ones(Nr+2,1);
  CellSize.y= [0.0];
  CellSize.z= [0.0];
  CellLocation.x= [1:Nr]'*dx-dx/2;
  CellLocation.y= [0.0];
  CellLocation.z= [0.0];
  FaceLocation.x= [0:Nr]'*dx;
  FaceLocation.y= [0.0];
  FaceLocation.z= [0.0];
elseif nargin==1
  % nonuniform 1D mesh
  facelocationR=varargin{1};
  n=size(facelocationR);
  if n(1)==1
      facelocationR=facelocationR';
  end
  Nr = length(facelocationR)-1;
  CellSize.x= [facelocationR(2)-facelocationR(1); ...
    facelocationR(2:end)-facelocationR(1:end-1); ...
    facelocationR(end)-facelocationR(end-1)];
  CellSize.y= [0.0];
  CellSize.z= [0.0];
  CellLocation.x= 0.5*(facelocationR(2:end)+facelocationR(1:end-1));
  CellLocation.y= [0.0];
  CellLocation.z= [0.0];
  FaceLocation.x= facelocationR;
  FaceLocation.y= [0.0];
  FaceLocation.z= [0.0];
end

MS=MeshCylindrical1D( ...
  [Nr,1], ...
  CellSize, ...
  CellLocation, ...
  FaceLocation, ...
  [1], ...
  [1]);
