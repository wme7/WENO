function outputf = WENO7_reconstruction(a,qo,ind)
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
%       |______________S0_______________|
%       |                               |
%       |       |______________S1_______________|
%       |       |                               |
%       |       |       |______________S2_______________|
%       |       |       |                               |
%       |       |       |       |_______________S3______________|
%     ..|---o---|---o---|---o---|---o---|---o---|---o---|---o---|...
%       | I{i-3}| I{i-2}| I{i-1}|  I{i} | I{i+1}| I{i+2}|| I{i+3}
%                               |+
%                             i-1/2
%
% WENO stencil: S{i} = [ I{i-3},...,I{i+3} ]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note: by using circshift over our domain, we are implicitly creating
% faqorable code that includes peridical boundary conditions. MD. 12.2012.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set Variables
qmmm=circshift(qo,ind*3);
qmm =circshift(qo,ind*2);
qm  =circshift(qo,ind);
qp  =circshift(qo,-ind);
qpp =circshift(qo,-ind*2);
qppp=circshift(qo,-ind*3);

% Reconstruction Polynomials
p0 = (-3*qmmm + 13*qmm - 23*qm  + 25*qo  )/12;
p1 = ( 1*qmm  - 5*qm   + 13*qo  +  3*qp  )/12;
p2 = (-1*qm   + 7*qo   +  7*qp  -  1*qpp )/12;
p3 = ( 3*qo   + 13*qp  -  5*qpp +  1*qppp)/12;

% Smooth Indicators
B0 = qm.*(134241*qm-114894*qo)   +qmmm.*(56694*qm-47214*qmm+6649*qmmm-22778*qo)...
    +25729*qo.^2  +qmm.*(-210282*qm+85641*qmm+86214*qo);
B1 = qo.*(41001*qo-30414*qp)     +qmm.*(-19374*qm+3169*qmm+19014*qo-5978*qp)...
    +6649*qp.^2   +qm.*(33441*qm-70602*qo+23094*qp);
B2 = qp.*(33441*qp-19374*qpp)    +qm.*(6649*qm-30414*qo+23094*qp-5978*qpp)...
    +3169*qpp.^2  +qo.*(41001*qo-70602*qp+19014*qpp);
B3 = qpp.*(85641*qpp-47214*qppp) +qo.*(25729*qo-114894*qp+86214*qpp-22778*qppp)...
    +6649*qppp.^2 +qp.*(134241*qp-210282*qpp+56694*qppp);

% Constants
epsilon = 1e-6;
g0 = 1/20; g1 = 9/20; g2 = 9/20; g3 = 1/20;

% Alpha weights
alpha0 = g0./(epsilon + B0/2880).^2;
alpha1 = g1./(epsilon + B1/2880).^2;
alpha2 = g2./(epsilon + B2/2880).^2;
alpha3 = g3./(epsilon + B3/2880).^2;
alphasum = alpha0 + alpha1 + alpha2 + alpha3;

% Non-linear weigths
w0 = alpha0./alphasum;
w1 = alpha1./alphasum;
w2 = alpha2./alphasum;
w3 = alpha3./alphasum;

% Numerical Flux at cell boundary, $u_{i+1/2}^{-}$;
qp = w0.*p0 + w1.*p1 + w2.*p2 + w3.*p3;

% Compute flux
fp = qp.*a; outputf = -(fp-circshift(fp,ind));
% End of reconstruction.
