close all;
clear variables;

constants; % Loads constants

%% Problem definition.
x = (0:0.05:10)';
finalTime = 1000e-9;

% Materials.

% Boundary conditions.

% Sources.
excPoint = floor(length(x)/2);
delay = 8e-9;
spread = 2e-9;

% Output requests.

% Initial fields.
initialEz = exp(- (x-5).^2);

%% Inits spatial semi-discretization.
cells = size(x,1);

ez=zeros(size(x,1),2);
hy=zeros(size(x,1),2);
if (exist('initialEz','var'))
    ez(:,1) = initialEz(:);
end

% Inits Absorbing Boundary Conditions.
exm2=0;             
exm1=0;
exn2=0;
exn1=0;

% Determines recursion coefficients
cfl = 1;
dx = sum(x(2:end)-x(1:(end-1)))/(length(x)-1);
dt = cfl*dx/c0;

cE = dt/eps0/dx;
cH = dt/mu0/dx;

%% Performs time integration.
tic;
for t=0:dt:finalTime
    % --- Updates E field ---
    for i=2:cells
        ez(i,2)=ez(i,1)+cE*(hy(i-1,1)-hy(i,1));
    end
%     ez(2:end,2) = ez(2:end,1) + cE.* (hy(1:(end-1),1)-hy(2:end,1));
    
    % --- Sources ---
%     ez(excPoint,2) = ez(excPoint,2) + exp(- 0.5*((t-delay)/spread)^2);
    
    % --- Boundary conditions ---
    ez(    1, 2) = -ez(      2, 2);
    ez(cells, 2) = -ez(cells-1, 2); % PEC
    
    % --- Updates H field ---
    for i=1:cells-1
        hy(i,2)=hy(i,1)+cH*(ez(i,2)-ez(i+1,2));
    end
%     hy(1:(end-1),2) = hy(1:(end-1),2) + cH.* (ez(1:(end-1),2) - ez(2:end,2));
    
    ez(:,1)=ez(:,2);
    hy(:,1)=hy(:,2);
    
    % --- Output requests ---
    subplot(2,1,1);
    hold off;
    plot(ez(:,2));
    hold on;
    axis([0 cells -1 1]);
    title(sprintf('FDTD Time = %.2f nsec',t*1e9))
    subplot(2,1,2);
    hold off;
    plot(hy(:,2));
    hold on;
    axis([0 cells -0.005 0.005]);
    
    drawnow;
end
toc;

