function [BCMatrix, BCRHS] = boundaryConditionSpherical3D(BC)
% It creates the matrix of coefficients based on the BC structure provided
% by the user for a spherical (r, theta, phi) 3D mesh. It also generates
% the right hand side vector of the linear system of equations.
% Here x=r, y=theta (polar angle), z=phi (azimuthal angle). The r and
% theta boundaries are identical to the cylindrical 3D case; the phi
% (front/back) boundary flux picks up the extra 1/(r*sin(theta)) metric
% factor of the spherical gradient.
%
% SYNOPSIS:
%       [BCMatrix, BCRHS] = boundaryConditionSpherical3D(BC)
%
% SEE ALSO: boundaryConditionCylindrical3D, boundaryCondition

% extract data from the mesh structure
Nxyz = BC.domain.dims;
Nx = Nxyz(1); Ntetta = Nxyz(2); Nz = Nxyz(3);
G=reshape(1:(Nx+2)*(Ntetta+2)*(Nz+2), Nx+2, Ntetta+2, Nz+2);
dx_1 = BC.domain.cellsize.x(1);
dx_end = BC.domain.cellsize.x(end);
dtetta_1 = BC.domain.cellsize.y(1);
dtetta_end = BC.domain.cellsize.y(end);
dz_1 = BC.domain.cellsize.z(1);
dz_end = BC.domain.cellsize.z(end);
rp = repmat(BC.domain.cellcenters.x, 1, Nz);
% metric factor r*sin(theta) evaluated on the (r, theta) faces used by the
% phi (front/back) boundaries
rp_rt = repmat(BC.domain.cellcenters.x, 1, Ntetta);
thetap_rt = repmat(BC.domain.cellcenters.y', Nx, 1);
sphf = rp_rt.*sin(thetap_rt);

% number of boundary nodes (exact number is 2[(m+1)(n+1)*(n+1)*(p+1)+(m+1)*p+1]:
nb = 8*((Nx+1)*(Ntetta+1)+(Nx+1)*(Nz+1)+(Ntetta+1)*(Nz+1));

% define the vectors to be used for the creation of the sparse matrix
ii = zeros(nb,1);
jj = zeros(nb,1);
s = zeros(nb,1);

% define the RHS column vector
BCRHS = zeros((Nx+2)*(Ntetta+2)*(Nz+2), 1);

% assign value to the corner nodes (useless cells)
q = 1:8;
ii(q) = BC.domain.corners; jj(q) = BC.domain.corners;
s(q) = 1; BCRHS(BC.domain.corners) = 0;

% assign values to the edges (useless cells)
q = q(end)+(1:length(BC.domain.edges));
ii(q) = BC.domain.edges; jj(q) = BC.domain.edges;
s(q) = 1; BCRHS(BC.domain.edges) = 0;

% Assign values to the boundary condition matrix and the RHS vector based
% on the BC structure
if (BC.top.periodic ==0) && (BC.bottom.periodic == 0)
    % top boundary
    j=Ntetta+2;
    i=2:Nx+1;
    k=2:Nz+1;
    q = q(end)+(1:Nx*Nz);
    ii(q) = G(i,j,k);  jj(q) = G(i,j,k);  s(q) = BC.top.b/2 + BC.top.a./(dtetta_end*rp);
    q = q(end)+(1:Nx*Nz);
    ii(q) = G(i,j,k);  jj(q) = G(i,j-1,k); s(q) = BC.top.b/2 - BC.top.a./(dtetta_end*rp);
    BCRHS(G(i,j,k)) = BC.top.c;

    % Bottom boundary
    j=1;
    i=2:Nx+1;
    k=2:Nz+1;
    q = q(end)+(1:Nx*Nz);
    ii(q) = G(i,j,k);  jj(q) = G(i,j+1,k);  s(q) = -(BC.bottom.b/2 + BC.bottom.a./(dtetta_1*rp)); % consider the reverse direction of normal
    q = q(end)+(1:Nx*Nz);
    ii(q) = G(i,j,k);  jj(q) = G(i,j,k); s(q) = -(BC.bottom.b/2 - BC.bottom.a./(dtetta_1*rp)); % consider the reverse direction of normal
    BCRHS(G(i,j,k)) = -(BC.bottom.c);
elseif (BC.top.periodic ==1) || (BC.bottom.periodic == 1) % periodic
    % top boundary
    j=Ntetta+2;
    i=2:Nx+1;
    k=2:Nz+1;
    q = q(end)+(1:Nx*Nz);
    ii(q) = G(i,j,k);  jj(q) = G(i,j,k);  s(q) = 1;
    q = q(end)+(1:Nx*Nz);
    ii(q) = G(i,j,k);  jj(q) = G(i,j-1,k);  s(q) = -1;
    q = q(end)+(1:Nx*Nz);
    ii(q) = G(i,j,k);  jj(q) = G(i,1,k); s(q) = dtetta_end/dtetta_1;
    q = q(end)+(1:Nx*Nz);
    ii(q) = G(i,j,k);  jj(q) = G(i,2,k); s(q) = -dtetta_end/dtetta_1;
    BCRHS(G(i,j,k)) = 0;

    % Bottom boundary
    j=1;
    i=2:Nx+1;
    k=2:Nz+1;
    q = q(end)+(1:Nx*Nz);
    ii(q) = G(i,j,k);  jj(q) = G(i,j,k);  s(q) = 1;
    q = q(end)+(1:Nx*Nz);
    ii(q) = G(i,j,k);  jj(q) = G(i,j+1,k);  s(q) = 1;
    q = q(end)+(1:Nx*Nz);
    ii(q) = G(i,j,k);  jj(q) = G(i,Ntetta+1,k); s(q) = -1;
    q = q(end)+(1:Nx*Nz);
    ii(q) = G(i,j,k);  jj(q) = G(i,Ntetta+2,k); s(q) = -1;
    BCRHS(G(i,j,k)) = 0;
end

if (BC.right.periodic == 0) && (BC.left.periodic == 0)
    % Right boundary
    i=Nx+2;
    j=2:Ntetta+1;
    k=2:Nz+1;
    q = q(end)+(1:Ntetta*Nz);
    ii(q) = G(i,j,k);  jj(q) = G(i,j,k);  s(q) = BC.right.b/2 + BC.right.a/dx_end;
    q = q(end)+(1:Ntetta*Nz);
    ii(q) = G(i,j,k);  jj(q) = G(i-1,j,k); s(q) = BC.right.b/2 - BC.right.a/dx_end;
    BCRHS(G(i,j,k)) = BC.right.c;

    % Left boundary
    i = 1;
    j=2:Ntetta+1;
    k=2:Nz+1;
    q = q(end)+(1:Ntetta*Nz);
    ii(q) = G(i,j,k);  jj(q) = G(i+1,j,k);  s(q) = -(BC.left.b/2 + BC.left.a/dx_1); % consider the reverse direction of normal
    q = q(end)+(1:Ntetta*Nz);
    ii(q) = G(i,j,k);  jj(q) = G(i,j,k); s(q) = -(BC.left.b/2 - BC.left.a/dx_1); % consider the reverse direction of normal
    BCRHS(G(i,j,k)) = -(BC.left.c);
elseif (BC.right.periodic == 1) || (BC.left.periodic == 1) % periodic
    % Right boundary
    i=Nx+2;
    j=2:Ntetta+1;
    k=2:Nz+1;
    q = q(end)+(1:Ntetta*Nz);
    ii(q) = G(i,j,k);  jj(q) = G(i,j,k);  s(q) = 1;
    q = q(end)+(1:Ntetta*Nz);
    ii(q) = G(i,j,k);  jj(q) = G(i-1,j,k);  s(q) = -1;
    q = q(end)+(1:Ntetta*Nz);
    ii(q) = G(i,j,k);  jj(q) = G(1,j,k); s(q) = dx_end/dx_1;
    q = q(end)+(1:Ntetta*Nz);
    ii(q) = G(i,j,k);  jj(q) = G(2,j,k); s(q) = -dx_end/dx_1;
    BCRHS(G(i,j,k)) = 0;

    % Left boundary
    i = 1;
    j=2:Ntetta+1;
    k=2:Nz+1;
    q = q(end)+(1:Ntetta*Nz);
    ii(q) = G(i,j,k);  jj(q) = G(i,j,k);  s(q) = 1;
    q = q(end)+(1:Ntetta*Nz);
    ii(q) = G(i,j,k);  jj(q) = G(i+1,j,k);  s(q) = 1;
    q = q(end)+(1:Ntetta*Nz);
    ii(q) = G(i,j,k);  jj(q) = G(Nx+1,j,k); s(q) = -1;
    q = q(end)+(1:Ntetta*Nz);
    ii(q) = G(i,j,k);  jj(q) = G(Nx+2,j,k); s(q) = -1;
    BCRHS(G(i,j,k)) = 0;
end

if (BC.front.periodic == 0) && (BC.back.periodic == 0)
    % Front boundary (phi+): flux carries the 1/(r*sin(theta)) metric
    k=Nz+2;
    i = 2:Nx+1;
    j=2:Ntetta+1;
    q = q(end)+(1:Nx*Ntetta);
    ii(q) = G(i,j,k);  jj(q) = G(i,j,k);  s(q) = BC.front.b/2 + BC.front.a./(sphf*dz_end);
    q = q(end)+(1:Nx*Ntetta);
    ii(q) = G(i,j,k);  jj(q) = G(i,j,k-1); s(q) = BC.front.b/2 - BC.front.a./(sphf*dz_end);
    BCRHS(G(i,j,k)) = BC.front.c;

    % Back boundary (phi-)
    k=1;
    i = 2:Nx+1;
    j=2:Ntetta+1;
    q = q(end)+(1:Nx*Ntetta);
    ii(q) = G(i,j,k);  jj(q) = G(i,j,k+1);  s(q) = -(BC.back.b/2 + BC.back.a./(sphf*dz_1)); % consider the reverse direction of normal
    q = q(end)+(1:Nx*Ntetta);
    ii(q) = G(i,j,k);  jj(q) = G(i,j,k); s(q) = -(BC.back.b/2 - BC.back.a./(sphf*dz_1)); % consider the reverse direction of normal
    BCRHS(G(i,j,k)) = -(BC.back.c);
elseif (BC.front.periodic == 1) || (BC.back.periodic == 1) % periodic
    % Front boundary
    k=Nz+2;
    i = 2:Nx+1;
    j=2:Ntetta+1;
    q = q(end)+(1:Nx*Ntetta);
    ii(q) = G(i,j,k);  jj(q) = G(i,j,k);  s(q) = 1;
    q = q(end)+(1:Nx*Ntetta);
    ii(q) = G(i,j,k);  jj(q) = G(i,j,k-1);  s(q) = -1;
    q = q(end)+(1:Nx*Ntetta);
    ii(q) = G(i,j,k);  jj(q) = G(i,j,1); s(q) = dz_end/dz_1;
    q = q(end)+(1:Nx*Ntetta);
    ii(q) = G(i,j,k);  jj(q) = G(i,j,2); s(q) = -dz_end/dz_1;
    BCRHS(G(i,j,k)) = 0;

    % Back boundary
    k=1;
    i = 2:Nx+1;
    j=2:Ntetta+1;
    q = q(end)+(1:Nx*Ntetta);
    ii(q) = G(i,j,k);  jj(q) = G(i,j,k);  s(q) = 1;
    q = q(end)+(1:Nx*Ntetta);
    ii(q) = G(i,j,k);  jj(q) = G(i,j,k+1);  s(q) = 1;
    q = q(end)+(1:Nx*Ntetta);
    ii(q) = G(i,j,k);  jj(q) = G(i,j,Nz+1); s(q) = -1;
    q = q(end)+(1:Nx*Ntetta);
    ii(q) = G(i,j,k);  jj(q) = G(i,j,Nz+2); s(q) = -1;
    BCRHS(G(i,j,k)) = 0;
end

% Build the sparse matrix of the boundary conditions
BCMatrix = sparse(ii(1:q(end)), jj(1:q(end)), s(1:q(end)), ...
    (Nx+2)*(Ntetta+2)*(Nz+2), (Nx+2)*(Ntetta+2)*(Nz+2));
end
