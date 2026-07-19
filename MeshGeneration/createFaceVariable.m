function facevar = createFaceVariable(meshvar, faceval)
% this function creates a face variable based on the geometry and mesh
% size. For instance it can be used to create a gravity field.
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
d = meshvar.dimension;
mn = meshvar.dims;

if (d ==1)
	xvalue = faceval(1).*ones(mn(1)+1,1);
    yvalue=[];
    zvalue=[];
elseif (d == 2)
    if numel(faceval)==2
        xvalue = faceval(1).*ones(mn(1)+1, mn(2));
        yvalue = faceval(2).*ones(mn(1), mn(2)+1);
        zvalue=[];
    else
        xvalue = faceval(1).*ones(mn(1)+1, mn(2));
        yvalue = faceval(1).*ones(mn(1), mn(2)+1);
        zvalue=[];
    end
elseif (d == 3)
    if numel(faceval)==3
        xvalue = faceval(1).*ones(mn(1)+1, mn(2), mn(3));
        yvalue = faceval(2).*ones(mn(1), mn(2)+1, mn(3));
        zvalue = faceval(3).*ones(mn(1), mn(2), mn(3)+1);
    else
        xvalue = faceval(1).*ones(mn(1)+1, mn(2), mn(3));
        yvalue = faceval(1).*ones(mn(1), mn(2)+1, mn(3));
        zvalue = faceval(1).*ones(mn(1), mn(2), mn(3)+1);
    end
end
facevar=FaceVariable(meshvar, xvalue, yvalue, zvalue);
