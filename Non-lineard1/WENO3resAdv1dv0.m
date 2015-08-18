function res = WENO3resAdv1dv0(u,flux,dflux,S,dx)
% Compute nvmerical fluxes at cell 'i' interfaces.
% Input:  u: cells average values,
%         flux: f(u),
%         dflux: df(u)/dx,
%         S: S(u), source term,
%         dx: cell size.
% Output: res = df/dx.
%
% coded by Manuel Diaz, NTU, 2013.10.20
%
% Domain cells (I{i}) reference:
%
%                |           |   u(i)    |           |
%                |  u(i-1)   |___________|           |
%                |___________|           |   u(i+1)  |
%                |           |           |___________|
%             ...|-----0-----|-----0-----|-----0-----|...
%                |    i-1    |     i     |    i+1    |
%                |-         +|-         +|-         +|
%              i-3/2       i-1/2       i+1/2       i+3/2
%
% ENO stencils (S{r}) reference:
%
%
%                               |___________S2__________|
%                               |                       |
%                       |___________S1__________|       |
%                       |                       |       |
%               |___________S0__________|       |       |
%             ..|---o---|---o---|---o---|---o---|---o---|...
%               | I{i-2}| I{i-1}|  I{i} | I{i+1}| I{i+2}|
%                                      -|
%                                     i+1/2
%
%
%               |___________S0__________|
%               |                       |
%               |       |___________S1__________|
%               |       |                       |
%               |       |       |___________S2__________|
%             ..|---o---|---o---|---o---|---o---|---o---|...
%               | I{i-2}| I{i-1}|  I{i} | I{i+1}| I{i+2}|
%                               |+
%                             i-1/2
%
% WENO stencil: S{i} = [ I{i-2},...,I{i+2} ]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% based on: http://www.mathworks.com/matlabcentral/fileexchange/40956
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N=numel(u);	i=2:(N+1);

% Preallocate memory and add ghost nodes at the ends
U(1) = u(1); U(i) = u(i-1); U(N+2) = u(end); U(N+3) = u(end);
alpha1 = zeros(size(U)); alpha2 = zeros(size(U)); 
R_minus = zeros(size(U)); R_plus = zeros(size(U)); 

% Lax-Friedrichs (LF) flux splitting
a=max(abs(dflux(U))); Fp = 0.5*(flux(U)+a*U); Fm = 0.5*(flux(U)-a*U);

% WENO3 "right" flux reconstruction
alpha1(i) = (1/3)./(eps + (Fp(i)-Fp(i-1)).^2).^2;   alpha2(i) = (2/3)./(eps + (Fp(i+1)-Fp(i)).^2).^2;
omega1 = alpha1./(alpha1 + alpha2);  omega2 = alpha2./(alpha1 + alpha2);
R_plus(i) = omega1(i).*(3/2*Fp(i) - 1/2*Fp(i-1)) + omega2(i).*(1/2*Fp(i) + 1/2*Fp(i+1));

% WENO3 "left" flux reconstruction
alpha1(i) = (1/3)./(eps + (Fm(i+2)-Fm(i+1)).^2).^2;   alpha2(i) = (2/3)./(eps + (Fm(i+1)-Fm(i)).^2).^2;   
omega1 = alpha1./(alpha1 + alpha2);  omega2 = alpha2./(alpha1 + alpha2);
R_minus(i) = omega1(i).*(3/2*Fm(i+1) - 1/2*Fm(i+2)) + omega2(i).*(1/2*Fm(i) + 1/2*Fm(i+1));

% Compute Residual
res(i-1) = (R_plus(i)+R_minus(i)-R_plus(i-1)-R_minus(i-1))/dx - S(u);