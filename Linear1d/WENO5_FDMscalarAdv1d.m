%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          Solving 1-D wave equation with Methods of lines and 
%       3rd order Weighted Essentially Non-Oscilaroty (MOL-WENO5)
%
%                du/dt + df/dx = S,  for x \in [a,b]
%                  where f = f(u): linear/nonlinear
%                     and S = s(u): source term
%
%             coded by Manuel Diaz, manuel.ade'at'gmail.com 
%              Institute of Applied Mechanics, 2012.08.20
%                               
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ref: G. Jiang and C.-W. Shu, Efficient implementation of weighted ENO
% schemes, Journal of Computational Physics, 126:202-228, (1996). 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Notes: Finite Difference implementation of WENO5 with SSP-RK33 time
% integration method. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; %close all; clc;

%% Parameters
    c = 1.0;   % advection velocity
   nx = 200;	% number of cells
  CFL = 0.4;    % Courant Number
 tEnd = 2.0;    % End time

%% Preprocess
      
% Build discrete domain
a=-1; b=1; dx=(b-a)/nx; x=a+dx/2:dx:b; 

% Build IC
ICcase=1;  % {1}Testing, {2}Costum ICs
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
t=0; it=0; u=u0; dt=dt0;

while t < tEnd
    % Update/correct time step
    if t+dt>tEnd, dt=tEnd-t; end
    
    % RK initialization
    uo = u;

    % 1st stage
    dF = resWENO5(u,c,dx);
    u = uo-dt*dF;
    
    % 2nd Stage
    dF = resWENO5(u,c,dx);
    u = 0.75*uo+0.25*(u-dt*dF);

    % 3rd stage
    dF = resWENO5(u,c,dx);
    u = (uo+2*(u-dt*dF))/3;
    
    % Update time and counter
    t=t+dt; it=it+1; 
    
    % Plot solution
    if rem(it,10) == 0
        plot(x,u0,'-x',x,u,'.'); axis(plotrange); drawnow;
    end
end
%% Final Plot
plot(x,u0,'-x',x,u,'.'); axis(plotrange)
title('WENO5, cell averages plot','interpreter','latex','FontSize',18);
xlabel('$\it{x}$','interpreter','latex','FontSize',14);
ylabel({'$\it{u(x)}$'},'interpreter','latex','FontSize',14);