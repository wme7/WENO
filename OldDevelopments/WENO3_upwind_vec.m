function h = WENO3_upwind_vec(v)
% *************************************************************************
% This subroutine 'upwinds' our information from left to right.
% This is done by assuming u(i) are the cell's averages flux information in
% our domain (see domain ref), the cell boundary fluxes are computed using
% WENO3 (k=2) reconstruction. 
%
% Input: u(i) = [u(i-2) u(i-1) u(i) u(i+1) u(i+2)];
% Output: h(i) = $u_{i+1/2}^{-}$
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
%                               |______S1_______|
%                               |               |
%                       |______S0_______|       |
%             ..|---o---|---o---|---o---|---o---|---o---|...
%               | I{i-2}| I{i-1}|  I{i} | I{i+1}| I{i+2}|
%                                      -|
%                                     i+1/2
%
% WENO stencil: S{i} = [ I{i-2},...,I{i+2} ]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note: by using circshift over our domain, we are implicitly creating
% favorable code that includes peridical boundary conditions. MD. 12.2012.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Right Flux
% Choose the positive fluxes, 'v', to compute the left cell boundary flux:
% $u_{i+1/2}^{-}$
vm  = circshift(v,[0 1]);
vp  = circshift(v,[0 -1]);

% Polynomials
p0n = (-vm + 3*v)/2;
p1n = ( v  + vp )/2;

% Smooth Indicators, Beta factors
B0n = (vm-v).^2; 
B1n = (v-vp).^2;

% Constants
d0n = 1/3; d1n = 2/3; epsilon = 1E-6;

% Alpha weights 
alpha0n = d0n./(epsilon + B0n).^2;
alpha1n = d1n./(epsilon + B1n).^2;
alphasumn = alpha0n + alpha1n;

% ENO stencils weigths
w0n = alpha0n./alphasumn;
w1n = alpha1n./alphasumn;

% Numerical Flux at cell boundary, $u_{i+1/2}^{-}$;
hn = w0n.*p0n + w1n.*p1n;

% Numerical Flux
% We assume advection speed positive (a > 0) therefore wind blows from left
% to right. We would use u_{i+1/2}^{-}, (un), for the numerical flux;
h = hn;