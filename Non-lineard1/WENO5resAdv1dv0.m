function res = WENO5resAdv1dv0(u,flux,dflux,S,dx)
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
N=numel(u);	i=3:(N+2);

% Preallocate memory and add ghost nodes at the ends
U(1) = u(1); U(2) = u(1); U(i) = u(i-2); U(N+3) = u(end); U(N+4) = u(end); U(N+5) = u(end);
alpha1 = zeros(size(U)); alpha2 = zeros(size(U)); alpha3 = zeros(size(U));
beta1 = zeros(size(U)); beta2 = zeros(size(U)); beta3 = zeros(size(U));
R_minus = zeros(size(U)); R_plus = zeros(size(U)); 

% Lax-Friedrichs (LF) flux splitting
a=max(abs(dflux(U))); Fp = 0.5*(flux(U)+a*U); Fm = 0.5*(flux(U)-a*U);

% WENO5 "right" flux reconstruction
beta1(i) = (13.0/12.0)*(Fp(i) - 2.0*Fp(i+1) + Fp(i+2)).^2 ...
         + (1.0/4.0)*(3.0*Fp(i) - 4.0*Fp(i+1) + Fp(i+2)).^2;
beta2(i) = (13.0/12.0)*(Fp(i-1) - 2.0*Fp(i) + Fp(i+1)).^2 ... 
         + (1.0/4.0)*(Fp(i-1) - Fp(i+1)).^2;
beta3(i) = (13.0/12.0)*(Fp(i-2) - 2.0*Fp(i-1) + Fp(i)).^2 ...
         + (1.0/4.0)*(Fp(i-2) - 4.0*Fp(i-1) + 3.0*Fp(i)).^2;

alpha1(i) = (3.0/10.0)./(eps + beta1(i)).^2;
alpha2(i) = (3.0/5.0)./(eps + beta2(i)).^2;
alpha3(i) = (1.0/10.0)./(eps + beta3(i)).^2;

omega1 = alpha1./(alpha1 + alpha2 + alpha3);
omega2 = alpha2./(alpha1 + alpha2 + alpha3);
omega3 = alpha3./(alpha1 + alpha2 + alpha3);

R_plus(i) = omega1(i).*(1.0/3.0*Fp(i) + 5.0/6.0*Fp(i+1) - 1.0/6.0*Fp(i+2)) ...
          + omega2(i).*(-1.0/6.0*Fp(i-1) + 5.0/6.0*Fp(i) + 1.0/3.0*Fp(i+1)) ...
          + omega3(i).*(1.0/3.0*Fp(i-2) - 7.0/6.0*Fp(i-1) + 11.0/6.0*Fp(i));

% WENO5 "left" flux reconstruction
beta1(i) = (13.0/12.0)*(Fm(i+1) - 2.0*Fm(i+2) + Fm(i+3)).^2 ...
         + (1.0/4.0)*(3.0*Fm(i+1) - 4.0*Fm(i+2) + Fm(i+3)).^2;
beta2(i) = (13.0/12.0)*(Fm(i) - 2.0*Fm(i+1) + Fm(i+2)).^2 ...
         + (1.0/4.0)*(Fm(i) - Fm(i+2)).^2;
beta3(i) = (13.0/12.0)*(Fm(i-1) - 2.0*Fm(i) + Fm(i+1)).^2 ...
         + (1.0/4.0)*(Fm(i-1) - 4.0*Fm(i) + 3.0*Fm(i+1)).^2;

alpha1 = (1.0/10.0)./(eps + beta1).^2;
alpha2 = (3.0/5.0)./(eps + beta2).^2;
alpha3 = (3.0/10.0)./(eps + beta3).^2;

omega1 = alpha1./(alpha1 + alpha2 + alpha3);
omega2 = alpha2./(alpha1 + alpha2 + alpha3);
omega3 = alpha3./(alpha1 + alpha2 + alpha3);

R_minus(i) = omega1(i).*(1.0/3.0*Fm(i+3) - 7.0/6.0*Fm(i+2) + 11.0/6.0*Fm(i+1)) ...
           + omega2(i).*(-1.0/6.0*Fm(i+2) + 5.0/6.0*Fm(i+1) + 1.0/3.0*Fm(i)) ...
           + omega3(i).*(1.0/3.0*Fm(i+1) + 5.0/6.0*Fm(i) - 1.0/6.0*Fm(i-1));

% Compute residual
res(i-2) = (R_plus(i)+R_minus(i)-R_plus(i-1)-R_minus(i-1))/dx - S(u);