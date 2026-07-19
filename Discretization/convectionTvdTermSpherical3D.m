function [M, RHS, Mx, My, Mz, RHSx, RHSy, RHSz] = ...
    convectionTvdTermSpherical3D(u, phi, FL)
% This function uses the TVD scheme to discretize a 3D convection term
% \grad . (u \phi) on a spherical (r, theta, phi) mesh. It returns the
% implicit upwind matrix M together with the explicit deferred TVD
% correction RHS (and their r/theta/phi components).
%
% SYNOPSIS:
%   [M, RHS, Mx, My, Mz, RHSx, RHSy, RHSz] = ...
%       convectionTvdTermSpherical3D(u, phi, FL)
%
% SEE ALSO: convectionTvdTermCylindrical3D, convectionUpwindTermSpherical3D,
%   convectionTvdRHSSpherical3D

% implicit upwind part
[M, Mx, My, Mz] = convectionUpwindTermSpherical3D(u);
% explicit TVD deferred-correction part
[RHS, RHSx, RHSy, RHSz] = convectionTvdRHSSpherical3D(u, phi, FL);
end
