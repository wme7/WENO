%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              Solving 1-D wave equation with 7th order
%          Weighted Essentially Non-Oscilaroty (MOL-WENO5-LF)
%
%                 du/dt + df/dx = S, for x \in [a,b]
%                  where f = f(u): linear/nonlinear
%                     and S = s(u): source term
%
%             coded by Manuel Diaz, manuel.ade'at'gmail.com 
%            Institute of Applied Mechanics, NTU, 2012.08.20
%                               
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ref: C.-W. Shu, High order weighted essentially non-oscillatory schemes
% for convection dominated problems, SIAM Review, 51:82-126, (2009). 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Notes: A fully conservative finite volume implementation of the method of
% lines (MOL) using WENO7 associated with SSP-RK45 time integration method. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Special Thx to Kang Wei-Yi for repoting bugs & corrections. M.D. 01.2015.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; %close all; clc;

%% Parameters
   nx = 0080;	% number of cells
  CFL = 0.40;	% Courant Number
 tEnd = 0.40;   % End time

fluxfun='buckley'; % select flux function
% Define our Flux function
switch fluxfun
    case 'linear'   % Scalar Advection, CFL_max: 0.65
        c=1; flux = @(w) c*w; 
        dflux = @(w) c*ones(size(w));
    case 'burgers' % Burgers, CFL_max: 0.40  
        flux = @(w) w.^2/2; 
        dflux = @(w) w; 
    case 'buckley' % Buckley-Leverett, CFL_max: 0.40 & tEnd: 0.40
        flux = @(w) 4*w.^2./(4*w.^2+(1-w).^2);
        dflux = @(w) 8*w.*(1-w)./(5*w.^2-2*w+1).^2;
end

sourcefun='dont'; % add source term
% Source term
switch sourcefun
    case 'add'
        S = @(w) 0.1*w.^2;
    case 'dont'
        S = @(w) zeros(size(w));
end

% Build discrete domain
a=-1; b=1; dx=(b-a)/nx; x=a+dx/2:dx:b; 

% Build IC
ICcase=2;  % {1}Testing, {2}Costum ICs
switch ICcase
    case 1 % Testing IC
        u0=TestingIC(x);  % Jiang and Shu IC
    case 2 % Guassian IC
        u0=CommonIC(x,9)-1; % cases 1-9 <- check them out!
    otherwise
        error('IC file not listed');
end

% Plot range
dl=0.1; plotrange=[a,b,min(u0)-dl,max(u0)+dl];

%% Solver Loop

% Low storage Runge-Kutta coefficients
rk4a = [            0.0 ...
        -567301805773.0/1357537059087.0 ...
        -2404267990393.0/2016746695238.0 ...
        -3550918686646.0/2091501179385.0  ...
        -1275806237668.0/842570457699.0];
rk4b = [ 1432997174477.0/9575080441755.0 ...
         5161836677717.0/13612068292357.0 ...
         1720146321549.0/2090206949498.0  ...
         3134564353537.0/4481467310338.0  ...
         2277821191437.0/14882151754819.0];
rk4c = [             0.0  ...
         1432997174477.0/9575080441755.0 ...
         2526269341429.0/6820363962896.0 ...
         2006345519317.0/3224310063776.0 ...
         2802321613138.0/2924317926251.0];

% Using a 4rd Order 5-stage SSPRK time integration
res_u = zeros(1,nx); % Runge-Kutta residual storage
     
% load initial conditions
t=0; it=0; u=u0;

while t < tEnd
	% Update/correct time step
    dt=CFL*dx/max(abs(u)); if t+dt>tEnd, dt=tEnd-t; end
    
	% Update time and iteration counter
    t=t+dt; it=it+1;
    
    % RK step
    for RKs = 1:5
        t_local = t + rk4c(RKs)*dt;
        dF = WENO7resAdv1d(u,flux,dflux,S,dx);
        res_u = rk4a(RKs)*res_u + dt*dF;
        u = u - rk4b(RKs)*res_u;
    end

    % Plot solution
    if rem(it,10) == 0
        plot(x,u0,'-x',x,u,'.'); axis(plotrange); 
        shg; drawnow;
    end
end

%% Final Plot
plot(x,u0,'-x',x,u,'.'); axis(plotrange);
title('WENO7, cell averages plot','interpreter','latex','FontSize',18);
xlabel('$\it{x}$','interpreter','latex','FontSize',14);
ylabel({'$\it{u(x)}$'},'interpreter','latex','FontSize',14);