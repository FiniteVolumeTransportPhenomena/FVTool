function [X, Y, Z] = cellLocations(m)
% X = cellLocations(m)
% this function returns the location of the cell centers as cell
% variables. It can later be used in defining properties that are variable
% in space.
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
%   faceLocations

% Written by Ali A. Eftekhari
% See the license file

N=m.dims;
X= createCellVariable(m, 0);
Y= createCellVariable(m, 0);
Z= createCellVariable(m, 0);
d = m.dimension;
switch d
case 1
    	X=createCellVariable(m, m.cellcenters.x);
    case 2
        X= createCellVariable(m, repmat(m.cellcenters.x, 1, N(2)));
        Y= createCellVariable(m, repmat(m.cellcenters.y', N(1), 1));
    case 3
        X= createCellVariable(m, repmat(m.cellcenters.x, 1, N(2), N(3)));
        Y= createCellVariable(m, repmat(m.cellcenters.y', N(1), 1, N(3)));
        z=zeros(1,1,N(3));
        z(1,1,:)= m.cellcenters.z;
        Z= createCellVariable(m, repmat(z, N(1), N(2), 1));
end
