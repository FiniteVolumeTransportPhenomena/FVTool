function cellnorm = normCellVector(cellvec)
% this function calculates the second norm of a cell vector field
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
d = cellvec.domain.dimension;

if (d ==1)
	cellnormval = abs(cellvec.xvalue);
elseif (d == 2)
	cellnormval = realsqrt(cellvec.xvalue.*cellvec.xvalue+ ...
        cellvec.yvalue.*cellvec.yvalue);
elseif (d == 3)
    cellnormval = realsqrt(cellvec.xvalue.*cellvec.xvalue+ ...
        cellvec.yvalue.*cellvec.yvalue+ ...
        cellvec.zvalue.*cellvec.zvalue);
end
BC = createBC(cellvec.domain);
c=cellBoundary(cellnormval, BC);
cellnorm=CellVariable(cellvec.domain, c);
