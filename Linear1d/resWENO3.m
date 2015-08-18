function dF = resWENO3(q,ax,dx)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Compute the Residual for 1d wave equation using WENO5
%
%                         Residual = dq/dx
%
%               coded by Manuel Diaz, NTU, 2012.08.20
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ax,zy: scalar advection velocities in x and y directions respectively.

% Along x
[axp,axm]=fluxsplitting_scalar(ax,'LF');
axm=circshift(axm,[0 1]);

dq=   WENO3_reconstruction(axm,q,[0 -1]);
dq=dq+WENO3_reconstruction(axp,q,[0 +1]);

% The Residual
dF = -dq/dx;