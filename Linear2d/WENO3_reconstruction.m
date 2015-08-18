function outputf = WENO3_reconstruction(a,qo,ind)
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
%                               |______S1_______|
%                               |               |
%                       |______S0_______|       |
%             ..|---o---|---o---|---o---|---o---|---o---|...
%               | I{i-2}| I{i-1}|  I{i} | I{i+1}| I{i+2}|
%                                      -|
%                                     i+1/2
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note: by using circshift over our domain, we are implicitly creating
% favorable code that includes peridical boundary conditions. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set Variables
qm =circshift(qo,ind);
qp =circshift(qo,-ind);

% Reconstruction Polynomials
up0 = (-qm + 3*qo)/2.;
up1 = ( qo +  qp )/2.;

% Smooth parameters
b0 = (qm-qo).^2; 
b1 = (qo-qp).^2;

% Constants
g0=1/3; g1=2/3; epsilon=1e-6;

% weigths
wt0 = g0./(epsilon+b0).^2;
wt1 = g1./(epsilon+b1).^2;
sum_wt=wt0+wt1;

% Non-linear weigths
w0 = wt0./sum_wt;
w1 = wt1./sum_wt;

% WENO polynomial
qp = w0.*up0 + w1.*up1;

% Compute flux
fp = qp.*a; outputf = -(fp-circshift(fp,ind));
% End of reconstruction.