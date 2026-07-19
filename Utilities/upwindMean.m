function phiFaceAverage = upwindMean(phi, u)
% This function gets the value of the field variable phi defined
% over the MeshStructure and calculates the upwind average on
% the cell faces, for a uniform mesh based on the direction of the velocity
% vector u
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

% extract data from the mesh structure

d = phi.domain.dimension;
if (d ==1)
	phiFaceAverage = upwindMean1D(phi, u);
elseif (d == 2)
	phiFaceAverage = upwindMean2D(phi, u);
elseif (d == 3)
    phiFaceAverage = upwindMean3D(phi, u);
end
