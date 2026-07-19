function [X, Y, Z] = faceLocations(m)
% [X, Y, Z] = faceLocations(m)
% m is a mesh type
% This function returns the X, Y, and Z
% each one is a face variable itself, and can be used for the calculation
% of face variable as a function of locations
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
%  cellLocations

% Written by Ali A. Eftekhari
% See the license file
X=createFaceVariable(m, 0);
Y=createFaceVariable(m, 0);
Z=createFaceVariable(m, 0);
N=m.dims;
d = m.dimension;
switch d
case 1
    	X.xvalue= m.facecenters.x;
    case 2
        X.xvalue= repmat(m.facecenters.x, 1, N(2));
        X.yvalue= repmat(m.cellcenters.y', N(1)+1, 1);
        Y.xvalue= repmat(m.cellcenters.x, 1, N(2)+1);
        Y.yvalue= repmat(m.facecenters.y', N(1), 1);
    case 3
        X.xvalue= repmat(m.facecenters.x, 1, N(2), N(3));
        X.yvalue= repmat(m.cellcenters.y', N(1)+1, 1, N(3));
        z=zeros(1,1,N(3));
        z(1,1,:)= m.cellcenters.z;
        X.zvalue= repmat(z, N(1)+1, N(2), 1);
        % Y
        Y.xvalue= repmat(m.cellcenters.x, 1, N(2)+1, N(3));
        Y.yvalue= repmat(m.facecenters.y', N(1), 1, N(3));
        Y.zvalue= repmat(z, N(1), N(2)+1, 1);
        % Z
        z=zeros(1,1,N(3)+1);
        z(1,1,:)= m.facecenters.z;
        Z.xvalue= repmat(m.cellcenters.x, 1, N(2), N(3)+1);
        Z.yvalue= repmat(m.cellcenters.y', N(1), 1, N(3)+1);
        Z.zvalue= repmat(z, N(1), N(2), 1);
end
