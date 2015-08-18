function dF = resWENO7(q,ax,ay,dx,strategy)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Compute the Residual for 2d wave equation using WENO5
%
%                    Residual: RHS of dq/dt
%
%              coded by Manuel Diaz, NTU, 2012.12.18
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Implementation based in ideas of adv2d.m by Guillaume Roullet, 2011.
% ref. http://stockage.univ-brest.fr/~roullet/codes.html
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ax,zy: scalar advection velocities in x and y directions respectively.

% Along x
[axp,axm]=fluxsplitting_scalar(ax,strategy);
axm=circshift(axm,[0 1]);

dq=   WENO7_reconstruction(axm,q,[0 -1]);
dq=dq+WENO7_reconstruction(axp,q,[0 +1]);

% Along y
[ayp,aym]=fluxsplitting_scalar(ay,strategy);
aym=circshift(aym,[1 0]);

dq=dq+WENO7_reconstruction(aym,q,[-1 0]);
dq=dq+WENO7_reconstruction(ayp,q,[+1 0]);

% The Residual
dF = -dq/dx;
