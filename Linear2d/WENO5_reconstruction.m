function outputf = WENO5_reconstruction(a,qo,ind)
% Compute nvmerical fluxes at cell 'i' interfaces.
% Input:  v(1:6) = positive fluxes - cells average values
%         u(1:6) = negative fluxes - cells average values
% Output: hn(1)  = Numerical flux @ x_{i+1/2}^{-} | right flux
%         hp(1)  = Numerical flux @ x_{i-1/2}^{+} | left flux
%
% coded by Manuel Diaz, NTU, 2012.12.20
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
% Note: by using circshift over our domain, we are implicitly creating
% favorable code that includes peridical boundary conditions. MD. 12.2012.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set Variables
qmm=circshift(qo,ind*2);
qm =circshift(qo,ind);
qp =circshift(qo,-ind);
qpp=circshift(qo,-ind*2);

% Reconstruction Polynomials
up0 = (2*qmm - 7*qm + 11*qo)/6.;
up1 = ( -qm  + 5*qo + 2*qp )/6.;
up2 = (2*qo  + 5*qp - qpp  )/6.;

% Smooth parameters
b0 = 13/12*(qmm-2*qm+qo ).^2 + 1/4*(qmm-4*qm+3*qo).^2; 
b1 = 13/12*(qm -2*qo+qp ).^2 + 1/4*(qm-qp).^2;
b2 = 13/12*(qo -2*qp+qpp).^2 + 1/4*(3*qo-4*qp+qpp).^2;

% Constants
g0=1/10; g1=6/10; g2=3/10; epsilon=1e-6;

% weigths
wt0 = g0./(epsilon+b0).^2;
wt1 = g1./(epsilon+b1).^2;
wt2 = g2./(epsilon+b2).^2;
sum_wt=wt0+wt1+wt2;

% Non-linear weigths
w0 = wt0./sum_wt;
w1 = wt1./sum_wt;
w2 = wt2./sum_wt;

% WENO5 polynomial
qp = w0.*up0 + w1.*up1 + w2.*up2;

% Compute flux
fp = qp.*a; outputf = -(fp-circshift(fp,ind));
% End of reconstruction.