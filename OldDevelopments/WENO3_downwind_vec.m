function h = WENO3_downwind_vec(u)
% *************************************************************************
% This subroutine 'downwinds' our information from right to left.
% This is done by assuming u(i) are the cell's averages flux information in
% our domain (see domain ref), the cell boundary fluxes are computed using
% WENO3 (k=2) reconstruction. 
%
% Input: u(i) = [u(i-2) u(i-1) u(i) u(i+1) u(i+2)];
% Output: h(i) = $u_{i-1/2}^{+}$
%
% Based on:
% C.W. Shu's Lectures notes on: 'ENO and WENO schemes for Hyperbolic
% Conservation Laws'  
%
% coded by Manuel Diaz, 02.10.2012, NTU Taiwan.
% *************************************************************************
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
%                       |______S1_______|
%                       |               |
%                       |       |______S0_______|
%             ..|---o---|---o---|---o---|---o---|---o---|...
%               | I{i-2}| I{i-1}|  I{i} | I{i+1}| I{i+2}|
%                               |+
%                             i-1/2
%
% WENO stencil: S{i} = [ I{i-2},...,I{i+2} ]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note: by using circshift over our domain, we are implicitly creating
% favorable code that includes peridical boundary conditions. MD. 12.2012.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Left Flux 
% Choose the negative fluxes, 'u', to compute the left cell boundary flux:
% $u_{i-1/2}^{+}$ 
um  = circshift(u,[0 1]);
up  = circshift(u,[0 -1]);

% Polynomials
p0p = ( um + u )/2;
p1p = (3*u - up)/2;

% Smooth Indicators, Beta factors
B0p = (um-u).^2; 
B1p = (u-up).^2;

% Constants
d0p = 2/3; d1p = 1/3; epsilon = 1E-6;

% Alpha weights 
alpha0p = d0p./(epsilon + B0p).^2;
alpha1p = d1p./(epsilon + B1p).^2;
alphasump = alpha0p + alpha1p;

% ENO stencils weigths
w0p = alpha0p./alphasump;
w1p = alpha1p./alphasump;

% Numerical Flux at cell boundary, $u_{i-1/2}^{+}$;
hp = w0p.*p0p + w1p.*p1p;

% Numerical Flux
% We assume negative advection speed,(a < 0),therefore wind blows from
% right to left. We would use u_{i-1/2}^{+}, (un), for the numerical flux;
h = hp;
