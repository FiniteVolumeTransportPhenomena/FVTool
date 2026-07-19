function [Mout, RHSout] = combineBC(BC, Meq, RHSeq)
%COMBINEBC This function combines the boundary condition equations with the
%main physical model equations, and delivers the matrix of coefficient and
%RHS to be solved for the internal cells. It is useful if one needs to use 
%and ODE solver for the accumulation term, i.e.
% d phi/ dt = M phi
%
% SYNOPSIS:
%    [Mout, RHSout] = combineBC(BC, Meq, RHSeq)
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


switch BC.domain.dimension
    case 1
        [Mout, RHSout] = combineBC1D(BC, Meq, RHSeq);
    case 2
        [Mout, RHSout] = combineBC2D(BC, Meq, RHSeq);
    case 3
        [Mout, RHSout] = combineBC3D(BC, Meq, RHSeq);
    otherwise
        error('FVTool:unsupportedGeometry', ...
            'combineBC: unsupported dimension %g', BC.domain.dimension);
end
