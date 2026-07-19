function RHSout = excludeGhostRHS(MS, RHS)
% this function cuts out the RHS values related to ghost cells and returns
% an RHS, which is only for internal cells
%
% SYNOPSIS:
%
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

% Written by Ali A. Eftekhari
% See the license file

% check the size of the variable and the mesh dimension
d = MS.dimension;

if (d ==1)
	Nx=MS.dims(1);
	G = 1:Nx+2;
	RHSout = RHS(reshape(G(2:end-1),Nx,1));
elseif (d == 2)
	Nxy = MS.dims;
	Nx = Nxy(1); Ny = Nxy(2);
	G=reshape(1:(Nx+2)*(Ny+2), Nx+2, Ny+2);
	RHSout = RHS(reshape(G(2:end-1,2:end-1),Nx*Ny,1));
elseif (d == 3)
	Nxyz = MS.dims;
	Nx = Nxyz(1); Ny = Nxyz(2); Nz = Nxyz(3);
	G=reshape(1:(Nx+2)*(Ny+2)*(Nz+2), Nx+2, Ny+2, Nz+2);
  RHSout = RHS(reshape(G(2:end-1,2:end-1,2:end-1),Nx*Ny*Nz,1));
end
