function phi_bc = BC2GhostCells(phi)
% This function assigns the value at the boundary to the ghost cells
%
% SYNOPSIS:
%   cell_var = BC2GhostCells(phi)
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

% Copyright (c) 2012-2016 Ali Akbar Eftekhari
% See the license file
d = phi.domain.dimension;
if (d ==1)
    phi.value(1) = 0.5*(phi.value(1)+phi.value(2));
    phi.value(end) = 0.5*(phi.value(end)+phi.value(end-1));
elseif (d == 2)
    phi.value(1,2:end-1) = 0.5*(phi.value(1,2:end-1)+phi.value(2,2:end-1));
    phi.value(end,2:end-1) = 0.5*(phi.value(end,2:end-1)+phi.value(end-1,2:end-1));
    phi.value(2:end-1, 1) = 0.5*(phi.value(2:end-1, 1)+phi.value(2:end-1, 2));
    phi.value(2:end-1, end) = 0.5*(phi.value(2:end-1, end)+phi.value(2:end-1, end-1));
elseif (d == 3)
    phi.value(1, 2:end-1, 2:end-1) = 0.5*(phi.value(1, 2:end-1, 2:end-1)+phi.value(2, 2:end-1, 2:end-1));
    phi.value(end, 2:end-1, 2:end-1) = 0.5*(phi.value(end, 2:end-1, 2:end-1)+phi.value(end-1, 2:end-1, 2:end-1));
    phi.value(2:end-1, 1, 2:end-1) = 0.5*(phi.value(2:end-1, 1, 2:end-1)+phi.value(2:end-1, 2, 2:end-1));
    phi.value(2:end-1, end, 2:end-1) = 0.5*(phi.value(2:end-1, end, 2:end-1)+phi.value(2:end-1, end-1, 2:end-1));
    phi.value(2:end-1, 2:end-1, 1) = 0.5*(phi.value(2:end-1, 2:end-1, 1)+phi.value(2:end-1, 2:end-1, 2));
    phi.value(2:end-1, 2:end-1, end) = 0.5*(phi.value(2:end-1, 2:end-1, end-1)+phi.value(2:end-1, 2:end-1, end));
end

phi_bc = phi;