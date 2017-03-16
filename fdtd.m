close all;
clear variables;

constants; % Loads constants

%% Problem definition.
L = 10;
x = (0:.025:L)';
finalTime = L/c0/4;

% Determines recursion coefficients
cfl = 1;
dx = sum(x(2:end)-x(1:(end-1)))/(length(x)-1);
dt = cfl*dx/c0;

cE = dt/eps0/dx;
cH = dt/mu0/dx;
% Materials.

% Boundary conditions.

% Sources.
% excPoint = floor(length(x)/4);
% scaPoint = floor(length(x)*3/4);
% delay = 8e-9;
% spread = 2e-9;

% Output requests.

% Initial fields.
spread = 0.5;
initialEz = analyticalGaussian(x,-dt/2,L,spread);


%% Inits spatial semi-discretization.
cells = size(x,1);

ez=zeros(size(x,1),2);
hy=zeros(size(x,1)-1,2);
if (exist('initialEz','var'))
    ez(:,1) = initialEz(:);
end
if (exist('initialHy','var'))
    hy(:,1) = initialHy(:);
end

% Inits Absorbing Boundary Conditions.
exm2=0;             
exm1=0;
exn2=0;
exn1=0;



% abc = zeros(2,1);

%% Performs time integration.
% tic;
for t=0:dt:(finalTime+dt/2)
    % --- Updates E field ---
    for i=2:(cells-1)
        ez(i,2)=ez(i,1)+cE*(hy(i-1,1)-hy(i,1));
    end
    
    % --- Sources ---
%     ez(excPoint,2) = ez(excPoint,2) + exp(- 0.5*((t-delay)/spread)^2);
%     
%     phaseShift = (x(scaPoint) - x(excPoint)) / c0;
%     ez(scaPoint,2) = ez(scaPoint,2) - ...
%         exp(- 0.5*((t-delay-phaseShift)/spread)^2);
    
    % --- Boundary conditions ---
    ez(    1, 2) = 0;
    ez(cells, 2) = 0; % PEC
%     ez(1,2)    =ez(    1,1) - 2*cE*hy(      1,1);
%     ez(cells,2)=ez(cells,1) + 2*cE*hy(cells-1,1); % PMC
%     ez(1,2) = ez(2,1) + ...
%         (c0*dt-dx)/(c0*dt+dx)*(ez(2,2) - ez(1,1));
%     ez(cells,2) = ez(cells-1,1) + ...
%         (c0*dt-dx)/(c0*dt+dx)*(ez(cells,2) - ez(cells-1,1)); % MUR
    
    % --- Updates H field ---
    for i=1:cells-1
        hy(i,2)=hy(i,1)+cH*(ez(i,2)-ez(i+1,2));
    end

%     hy(excPoint,2) = hy(excPoint,2) + ...
%         exp(- 0.5*((t+dt/2-delay)/spread)^2)/eta0;
%     hy(scaPoint,2) = hy(scaPoint,2) - ...
%         exp(- 0.5*((t+dt/2-delay-phaseShift)/spread)^2)/eta0;

    ez(:,1)=ez(:,2);
    hy(:,1)=hy(:,2);

%     % --- Output requests ---
%     subplot(2,1,1);
%     hold off;
%     plot(x,ez(:,1));
%     hold on;
%     xAn = 0:0.001:L;
%     plot(xAn,analyticalGaussian(xAn,t+dt/2,L,spread),'-r');
%     axis([x(1) x(end) -1 1]);
%     title(sprintf('FDTD Time = %.2f nsec',t*1e9))
%     subplot(2,1,2);
%     hold off;
%     plot(x(1:end-1),hy(:,1));
%     hold on;
%     axis([x(1) x(end) -0.005 0.005]);
%     pause(.0025);
%     drawnow;
end
% toc;

Eanalytical = analyticalGaussian(x,t+dt/2,L,spread);

fprintf('For dx=%e , L^2 error: %e\n', ...
    dx, sum(abs(ez(:,2)-Eanalytical))/cells);
