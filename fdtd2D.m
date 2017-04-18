close all;
clear variables;

constants; % Loads constants

%% Problem definition.
L = 10;
x = (0:.25:L)';
y = (0:.25:L)';
finalTime = L/c0;

% Determines recursion coefficients
cfl = 1;
dx = sum(x(2:end)-x(1:(end-1)))/(length(x)-1);
dy = sum(y(2:end)-y(1:(end-1)))/(length(y)-1);
%
dt = cfl/c0/sqrt((dx)^-2+(dy)^-2);

cEx = dt/eps0/dx;
cEy = dt/eps0/dy;

cHx = dt/mu0/dx;
cHy = dt/mu0/dy;

% Initial fields.
spread = 0.5;
initialEz = analyticalGaussian2D(x,y,-dt/2,L,spread);



%% Inits spatial semi-discretization.
cellsX = size(x,1);
cellsY = size(y,1);

ez=zeros(size(x,1),size(y,1),2);
hx=zeros(size(x,1)-1,size(y,1)-1,2);
hy=zeros(size(x,1)-1,size(y,1)-1,2);
if  (exist('initialEz','var'))
     ez(:,:,1) = initialEz(:,:);
end
if  (exist('initialHx','var'))
     hy(:,:,1) = initialHx(:,:);
end
if  (exist('initialHy','var'))
     hy(:,:,1) = initialHy(:,:);
end



xAn=x;
yAn=y;
[x,y]=meshgrid(x,y); % Mallado para las gráficas

%% Performs time integration.
% tic;
for t=0:dt:(finalTime+dt/2)
    % --- Updates E field ---
    for i=2:(cellsX-1)
        for j=2:(cellsY-1)
            ez(i,j,2)=ez(i,j,1)+cEx*(hy(i,j,1)-hy(i-1,j,1))-...
                      cEy*(hx(i,j,1)-hx(i,j-1,1));
        end
    end
    
    
    % --- Boundary conditions ---
    ez(1,1,2) = 0;
    ez(cellsX,1,2) = 0; %  condiciones PEC
    ez(1,cellsY,2) = 0;
    ez(cellsX,cellsY,2) = 0;

    
    % --- Updates H field ---
    for i=1:(cellsX-1)
        for j=1:(cellsY-1)
            hx(i,j,2)=hy(i,j,1)-cHy*(ez(i,j+1,2)-ez(i,j,2));
            hy(i,j,2)=hy(i,j,1)+cHx*(ez(i+1,j,2)-ez(i,j,2));
        end
    end
    
    ez(:,:,1)=ez(:,:,2);
    hx(:,:,1)=hx(:,:,2);
    hy(:,:,1)=hy(:,:,2);

% Representaciones gráficas

subplot(2,2,1)
hold off;
surf(x,y,ez(:,:,1))
title('Ez(x,y,t)')
hold on;

subplot(2,2,2)
hold off;
surf(x,y,analyticalGaussian2D(xAn,yAn,t+dt/2,L,spread));
title('Analyticalgaussian2D(x,y,t)')
hold on;

subplot(2,2,3)
hold off;
surf(x(1:end-1,1:end-1),y(1:end-1,1:end-1),hx(:,:,1))
title('Hx(x,y,t)')
hold on;

subplot(2,2,4)
hold off;
surf(x(1:end-1,1:end-1),y(1:end-1,1:end-1),hy(:,:,1))
title('Hy(x,y,t)')
hold on;

drawnow;
end
% toc;


