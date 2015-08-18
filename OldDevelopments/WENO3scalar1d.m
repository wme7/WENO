%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         Solving 1-D wave equation using Method of Lines and 
%       3rd order Weighted Essentially Non-Oscilaroty (MOL-WENO3)
%
%                du/dt + df/dx = 0,  for x \in [a,b]
%                 where f = c*u, and c is a scalar.
%
%             coded by Manuel Diaz, manuel.ade'at'gmail.com 
%              Institute of Applied Mechanics, 2012.08.20
%                               
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ref: C.-W. Shu, High order weighted essentially non-oscillatory schemes
%      for convection dominated problems, SIAM Review, 51:82-126, (2009). 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTE: A simple vectorized implementation without time integration scheme.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; %close all; clc;

%% Parameters
    c = -1.0;	% advection velocity
   nx = 0200;	% number of cells
  CFL = 0.30;	% Courant Number
 tEnd = 4.00;   % End time

%% Preprocess
      
% Build discrete domain
a=-1; b=1; dx=(b-a)/nx; x=a+dx/2:dx:b; 

% Build IC
ICcase=2;  % {1}Testing, {2}Costum ICs
switch ICcase
    case 1 % Testing IC
        u0=TestingIC(x);  % Jiang and Shu IC
    case 2 % Guassian IC
        u0=CommonIC(x,10)-1; % cases 1-9 <- check them out!
    otherwise
        error('IC file not listed');
end

% Time discretization
dt0=CFL*dx/abs(c);

% Plot range
dl=0.1; plotrange=[a,b,min(u0)-dl,max(u0)+dl];

%% Solver Loop

% load initial conditions
dt=dt0; t=0; it=0; u=u0;

while t < tEnd
	% Update/correct time step
    if t+dt>tEnd, dt=tEnd-t; end
    
    % Update time and iteration counter
    t=t+dt; it=it+1;
    
    % Initial step
    uo = u;
    
    if c > 0;
        % Upwinding WENO
        ur = WENO3_upwind_vec(u);
        h = [ur(end),ur]; % periodic BC
        uR = h(2:end); uL = h(1:end-1);
        
        % Compute solution of next time step using WENO Upwind
        u = uo - c*dt/dx*(uR-uL);

    else % c < 0;
        % Downwinding WENO 
        ul = WENO3_downwind_vec(u);
        h = [ul,ul(1)]; % periodic BC
        uR = h(2:end); uL = h(1:end-1);

        % Compute solution of next time step using WENO Upwind
        u = uo - c*dt/dx*(uR-uL);
    end
    
    % Plot solution
    if rem(it,10) == 0
        plot(x,u0,'-x',x,u,'.'); axis(plotrange); shg; drawnow;
    end
end
toc
%% Final Plot
plot(x,u0,'-x',x,u,'.'); axis(plotrange);
title('WENO3, cell averages plot','interpreter','latex','FontSize',18);
xlabel('$\it{x}$','interpreter','latex','FontSize',14);
ylabel({'$\it{u(x)}$'},'interpreter','latex','FontSize',14);