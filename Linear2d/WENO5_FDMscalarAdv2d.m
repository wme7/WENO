%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               Solving 2-D wave equation with 5th order
%             Weighted Essentially Non-Oscilaroty (WENO5)
%
%            dq/dt + df/dx + dg/dy = 0,  for x,y \in [a,b;c,d]
%                 where f = u*q  and  g = v*q
%
%              coded by Manuel Diaz, NTU, 2012.12.18
%                               
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Implementation based on ideas of adv2d.m by Guillaume Roullet, 2011.
% ref. http://stockage.univ-brest.fr/~roullet/codes.html
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; %close all; clc;

%% Parameters
     u = -0.5;  % Scalar velocity in x direction
     v =  1.0;  % Scalar velocity in y direction
   CFL =  0.6;  % CFL condition
  tEnd = 1.00;  % Final time

%% Preprocess
% Domain discretization
a=0; b=1; c=0; d=1; [x,dx,y,dy]=grid2d(a,b,80,c,d,80);

% set IC
q0=IC2d(x,y,2); %{1} 4 Quadrants, {2} Square Jump

% Time discretization
dt=min(dy,dx)*CFL/max(abs(v),abs(u)); t=0:dt:tEnd; 

% set plot range
plotrange=[a,b/dx,c,d/dx,min(min(q0)),max(max(q0))];
    
%% Solver Loop 
% load initial conditions
q=q0; it=0;

tic
for tstep=t
    % RK initialization
    qo=q;
       
    % 1st stage
    dF = resWENO5(q,u,v,dx,'LF');
    q = qo-dt*dF;
    
    % 2nd Stage
    dF = resWENO5(q,u,v,dx,'LF');
    q = 0.75*qo+0.25*(q-dt*dF);
    
    % 3rd stage
    dF = resWENO5(q,u,v,dx,'LF');
    q = (qo+2*(q-dt*dF))/3;
    
    % plot
    if(rem(it,10)==0)
        subplot(1,2,1); mesh(q); colormap Copper; axis(plotrange);
        %colorbar('location','EastOutside');
        title(['WENO5, dx = ',num2str(dx),', dy = ',num2str(dy),', time: ',num2str(tstep)])
        xlabel('x points'); ylabel('y points'); zlabel('q(x,y)');
        subplot(1,2,2); contourf(q); colormap Copper;
        title(['WENO5, dx = ',num2str(dx),', dy = ',num2str(dy),', time: ',num2str(tstep)])
        xlabel('x points'); ylabel('y points');
        drawnow
    end
end
toc