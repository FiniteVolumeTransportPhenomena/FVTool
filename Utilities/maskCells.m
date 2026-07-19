function [M_out, RHS_out] = maskCells(meshstruct, M, RHS, cellIndex, cellValue)
% This function masks the specified cells by giving them a constant value.
% It modifies the matrix of coefficient and the RHS vector
%
% SYNOPSIS:
%   [M_out, RHS_out] = maskCells(meshstruct, M, RHS, cellIndex, cellValue)
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

% extract some data
d = meshstruct.dimension;
domain_size = meshstruct.dims+2; % 2 is added for the ghost cells
M_size = size(M);
if d ==1
    i = sub2ind([domain_size 1], cellIndex(:,1));
elseif d==2
    i = sub2ind(domain_size, cellIndex(:,1), cellIndex(:,2));
elseif d==3
    i = sub2ind(domain_size, cellIndex(:,1), cellIndex(:,2), cellIndex(:,2));
end

% define the new masked matrix of coefficients
M_masked = sparse(i,i,1, M_size(1), M_size(2));
RHS_masked = zeros(length(RHS),1);
RHS_masked(i) = cellValue;

% zero the masked rows in the original matrix
M(i, :) = 0;
RHS(i) = 0;

% add the new masked matrix to the modifed original
M_out = M+M_masked;
RHS_out = RHS + RHS_masked;
