function dF = resWENO3(q,ax,ay,dx,strategy)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Compute the Residual for 2d wave equation using Upwind
%
%                    Residual: RHS of dq/dt
%
%              coded by Manuel Diaz, NTU, 2012.12.20
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ax,zy: scalar advection velocities in x and y directions respectively.

% Along x
[axp,axm]=fluxsplitting_scalar(ax,strategy);
axm=circshift(axm,[0 1]);

dq=   WENO3_reconstruction(axp,q,[0 -1]);
dq=dq+WENO3_reconstruction(axm,q,[0 +1]);

% Along y
[ayp,aym]=fluxsplitting_scalar(ay,strategy);
aym=circshift(aym,[1 0]);

dq=dq+WENO3_reconstruction(ayp,q,[-1 0]);
dq=dq+WENO3_reconstruction(aym,q,[+1 0]);

% The Residual
dF = -dq/dx;