function RHS = constantSourceTerm(phi)
% RHS vector for an explicit source term
% k is a cell variable
%
% SYNOPSIS:
%   RHS = constantSourceTerm(k)
%
% PARAMETERS:
%   phi: Cell Variable
%
% RETURNS:
%   RHS: vector
%
% EXAMPLE:
%
% SEE ALSO:
%

switch phi.domain.dimension
    case 1
        RHS = constantSourceTerm1D(phi);
    case 2
        RHS = constantSourceTerm2D(phi);
    case 3
        RHS = constantSourceTerm3D(phi);
    otherwise
        error('FVTool:unsupportedGeometry', ...
            'constantSourceTerm: unsupported dimension %g', phi.domain.dimension);
end
