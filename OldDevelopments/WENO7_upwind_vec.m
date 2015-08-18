function h = WENO7_upwind_vec(vo)
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
%                               |_______________S3______________|
%                               |                               |
%                       |______________S2_______________|       |
%                       |                               |       |
%               |______________S1_______________|       |       |
%               |                               |       |       |
%       |_______________S0______________|       |       |       |
%     ..|---o---|---o---|---o---|---o---|---o---|---o---|---o---|...
%       | I{i-3}| I{i-2}| I{i-1}|  I{i} | I{i+1}| I{i+2}| I{i+3}|
%                                      -|
%                                     i+1/2
%
% WENO stencil: S{i} = [ I{i-3},...,I{i+3} ]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note: by using circshift over our domain, we are implicitly creating
% favorable code that includes peridical boundary conditions. MD. 12.2012.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set Variables
vmmm=circshift(vo,[0 3]);
vmm =circshift(vo,[0 2]);
vm  =circshift(vo,[0 1]);
vp  =circshift(vo,[0 -1]);
vpp =circshift(vo,[0 -2]);
vppp=circshift(vo,[0 -3]);

% Reconstruction Polynomials
p0n = (-3*vmmm + 13*vmm - 23*vm  + 25*vo  )/12;
p1n = ( 1*vmm  - 5*vm   + 13*vo  +  3*vp  )/12;
p2n = (-1*vm   + 7*vo   +  7*vp  -  1*vpp )/12;
p3n = ( 3*vo   + 13*vp  -  5*vpp +  1*vppp)/12;

% Smooth Indicators
B0n = vm.*(134241*vm-114894*vo)   +vmmm.*(56694*vm-47214*vmm+6649*vmmm-22778*vo)...
        +25729*vo.^2  +vmm.*(-210282*vm+85641*vmm+86214*vo);
B1n = vo.*(41001*vo-30414*vp)     +vmm.*(-19374*vm+3169*vmm+19014*vo-5978*vp)...
        +6649*vp.^2   +vm.*(33441*vm-70602*vo+23094*vp);
B2n = vp.*(33441*vp-19374*vpp)    +vm.*(6649*vm-30414*vo+23094*vp-5978*vpp)...
        +3169*vpp.^2  +vo.*(41001*vo-70602*vp+19014*vpp);
B3n = vpp.*(85641*vpp-47214*vppp) +vo.*(25729*vo-114894*vp+86214*vpp-22778*vppp)...
        +6649*vppp.^2 +vp.*(134241*vp-210282*vpp+56694*vppp);

% Constants
epsilon = 1e-6;
g0 = 1/20; g1 = 9/20; g2 = 9/20; g3 = 1/20;

% Alpha weights
alpha0n = g0./(epsilon + B0n).^2;
alpha1n = g1./(epsilon + B1n).^2;
alpha2n = g2./(epsilon + B2n).^2;
alpha3n = g3./(epsilon + B3n).^2;
alphasvmn = alpha0n + alpha1n + alpha2n + alpha3n;

% Non-linear weigths
w0n = alpha0n./alphasvmn;
w1n = alpha1n./alphasvmn;
w2n = alpha2n./alphasvmn;
w3n = alpha3n./alphasvmn;

% Nvmerical Flux at cell boundary, $u_{i+1/2}^{-}$;
vn = w0n.*p0n + w1n.*p1n + w2n.*p2n + w3n.*p3n;

% Nvmerical Flux
% We assvme advection speed positive (a > 0) therefore wind blows from left
% to right. We would use u_{i+1/2}^{-}, (un), for the numerical flux;
h = vn;
