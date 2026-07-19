function faceGrad = gradientTermFixedBC(phi)
% this function calculates the gradient of parameter phi in x,y, and z 
% directions. It takes care of the often nonphysical
% values of the ghost cells. Note that phi is not a variable but a parameter calculated
% with a function over a domain. Make sure that phi is calculated by 
% BC2GhostCells (usually but not necessarily in combination with celleval); 
% otherwise, do not use this function as it leads to wrong 
% values at the boundaries.
% it checks for the availability of the ghost variables and use them, otherwise
% estimate them, assuming a zero gradient on the boundaries.
%
% SYNOPSIS:
%   faceGrad = gradientTermFixedBC(phi)
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
faceGrad = gradientTerm(phi);
d = phi.domain.dimension;
if (d == 1)
    faceGrad.xvalue(1) = 2*faceGrad.xvalue(1);
    faceGrad.xvalue(end) = 2*faceGrad.xvalue(end);
elseif (d == 2)
    faceGrad.xvalue(1,:) = 2*faceGrad.xvalue(1,:);
    faceGrad.xvalue(end,:) = 2*faceGrad.xvalue(end,:);
    faceGrad.yvalue(:,1) = 2*faceGrad.yvalue(:,1);
    faceGrad.yvalue(:,end) = 2*faceGrad.yvalue(:,end);
elseif (d == 3)
    faceGrad.xvalue(1,:,:) = 2*faceGrad.xvalue(1,:,:);
    faceGrad.xvalue(end,:,:) = 2*faceGrad.xvalue(end,:,:);
    faceGrad.yvalue(:,1,:) = 2*faceGrad.yvalue(:,1,:);
    faceGrad.yvalue(:,end,:) = 2*faceGrad.yvalue(:,end,:);
    faceGrad.zvalue(:,:,1) = 2*faceGrad.zvalue(:,:,1);
    faceGrad.zvalue(:,:,end) = 2*faceGrad.zvalue(:,:,end);
end
