function cellvec = createCellVector(meshvar, cellval)
% this function creates a vector field based on the geometry and mesh
% size. For instance it can be used to create a gravity field.
%
% SYNOPSIS:
%   cellvec = createCellVector(meshvar, cellval)
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
%     cellBoundary

% Written by Ali A. Eftekhari
% See the license file

% check the size of the variable and the mesh dimension
d = meshvar.dimension;
mn = meshvar.dims;

if (d ==1)
	xvalue = cellval(1).*ones(mn(1),1);
    yvalue=[];
    zvalue=[];
elseif (d == 2)
    if numel(cellval)==2
        xvalue = cellval(1).*ones(mn(1), mn(2));
        yvalue = cellval(2).*ones(mn(1), mn(2));
        zvalue=[];
    else
        xvalue = cellval(1).*ones(mn(1), mn(2));
        yvalue = cellval(1).*ones(mn(1), mn(2));
        zvalue=[];
    end
elseif (d == 3)
    if numel(cellval)==3
        xvalue = cellval(1).*ones(mn(1), mn(2), mn(3));
        yvalue = cellval(2).*ones(mn(1), mn(2), mn(3));
        zvalue = cellval(3).*ones(mn(1), mn(2), mn(3));
    else
        xvalue = cellval(1).*ones(mn(1), mn(2), mn(3));
        yvalue = cellval(1).*ones(mn(1), mn(2), mn(3));
        zvalue = cellval(1).*ones(mn(1), mn(2), mn(3));
    end
end
cellvec=CellVector(meshvar, xvalue, yvalue, zvalue);
