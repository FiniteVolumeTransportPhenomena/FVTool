function M = linearSourceTerm(k)
% Matrix of coefficients for a linear source term in the form of k \phi
%
% k is a cell variable
%
% SYNOPSIS:
%   M = linearSourceTerm(k)
%
% PARAMETERS:
%	k   - source term, CellVariable
%
% RETURNS:
%
%
% EXAMPLE:
%
% SEE ALSO:
%

switch k.domain.dimension
    case 1
        M = linearSourceTerm1D(k);
    case 2
        M = linearSourceTerm2D(k);
    case 3
        M = linearSourceTerm3D(k);
    otherwise
        error('FVTool:unsupportedGeometry', ...
            'linearSourceTerm: unsupported dimension %g', k.domain.dimension);
end
