function h = WENO7_downwind_vec(uo)
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
%       |______________S0_______________|
%       |                               |
%       |       |______________S1_______________|
%       |       |                               |
%       |       |       |______________S2_______________|
%       |       |       |                               |
%       |       |       |       |_______________S3______________|
%     ..|---o---|---o---|---o---|---o---|---o---|---o---|---o---|...
%       | I{i-3}| I{i-2}| I{i-1}|  I{i} | I{i+1}| I{i+2}| I{i+3}|
%                               |+
%                             i-1/2
%
% WENO stencil: S{i} = [ I{i-3},...,I{i+3} ]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note: by using circshift over our domain, we are implicitly creating
% favorable code that includes peridical boundary conditions. MD. 12.2012.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Set Variables
ummm=circshift(uo,[0 3]);
umm =circshift(uo,[0 2]);
um  =circshift(uo,[0 1]);
up  =circshift(uo,[0 -1]);
upp =circshift(uo,[0 -2]);
uppp=circshift(uo,[0 -3]);

%% Reconstruction Polynomials
p0p = ( 1*ummm - 5*umm  + 13*um  +  3*uo  )/12;
p1p = (-1*umm  + 7*um   +  7*uo  -  1*up  )/12;
p2p = ( 3*um   + 13*uo  -  5*up  +  1*upp )/12;
p3p = (25*uo   - 23*up  + 13*upp -  3*uppp)/12;

%% Smooth Indicators
B0p = um.*(134241*um-114894*uo)   +ummm.*(56694*um-47214*umm+6649*ummm-22778*uo)...
        +25729*uo.^2  +umm.*(-210282*um+85641*umm+86214*uo);
B1p = uo.*(41001*uo-30414*up)     +umm.*(-19374*um+3169*umm+19014*uo-5978*up)...
        +6649*up.^2   +um.*(33441*um-70602*uo+23094*up);
B2p = up.*(33441*up-19374*upp)    +um.*(6649*um-30414*uo+23094*up-5978*upp)...
        +3169*upp.^2  +uo.*(41001*uo-70602*up+19014*upp);
B3p = upp.*(85641*upp-47214*uppp) +uo.*(25729*uo-114894*up+86214*upp-22778*uppp)...
        +6649*uppp.^2 +up.*(134241*up-210282*upp+56694*uppp);

% Constants
epsilon = 1e-6;
g0 = 1/20; g1 = 9/20; g2 = 9/20; g3 = 1/20;

% Alpha weights
alpha0p = g0./(epsilon + B0p).^2;
alpha1p = g1./(epsilon + B1p).^2;
alpha2p = g2./(epsilon + B2p).^2;
alpha3p = g3./(epsilon + B3p).^2;
alphasump = alpha0p + alpha1p + alpha2p + alpha3p;

% Non-linear weigths
w0p = alpha0p./alphasump;
w1p = alpha1p./alphasump;
w2p = alpha2p./alphasump;
w3p = alpha3p./alphasump;

%% Numerical Flux at cell boundary, $u_{i+1/2}^{-}$;
up = w0p.*p0p + w1p.*p1p + w2p.*p2p + w3p.*p3p;

% Numerical Flux
% We assume negative advection speed,(a < 0),therefore wind blows from
% right to left. We would use u_{i-1/2}^{+}, (un), for the numerical flux;
h = up;
