close all;
clear variables;

constants; % Loads constants

%% Campo electromagnético + lámina delgada

% Este código tiene por objetivo estudiar lo mismo que hemos estudiado en
% clase, solo que ahora ponemos una lámina delgada. La idea es la misma,
% solo que hay que tener en cuenta que la la´mina tiene una conductividad
% eléctrica y un epsilon que no es el del vacío.

% Como podéis observar, el código es prácticamente el mismo (¿no os parece
% sospechoso las similitudes entre uno y otro?), la diferencia está en que
% hay que ubicar una lámina dieléctrica y conductora. Una vez ubicada,
% hemos de decirle al programa su esilon y su sigma.

% Esto se logra mediante unas M-funciones que he creado (algo mío tenía que
% haber en este trabajo) las cuales establecen la región en la que
% físicamente hay lámina y, fuera de ella, no hay sigma y la epsilon es la
% del vacío.

%% Problem definition.
x = (0:0.05:10)';
finalTime = 120e-9;

% Materials.

sigma1=0; % Lámina dieléctrica
% sigma1= 1.60*10^(-5); % Silicio
 sigma1= 3.78*10^(7); % Aluminio
epsilon1=10*eps0; % Silicio (supuestamente)
% epsilon1=100*eps0;
epsvacio=eps0;

% Boundary conditions.

% Sources.
excPoint = floor(length(x)/2);
scaPoint=floor(length(x)*3/4);
delay = 8e-9;
spread = 2e-9;

% Output requests.

% Initial fields.
initialEz = exp(- (x-1).^2/0.01);

%% Inits spatial semi-discretization.
cells = size(x,1);

ez=zeros(cells,2);
hy=zeros(cells,2);
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


%% Metiendo la lámina

%k1=floor((cells)/2);
%k2=cells;
k1=floor((cells-20)/2);
k2=cells+22;
k2=ceil(k2/2);

x1= x(k1); % Extremo izq. de la lámina
x2= x(k2); % Extremo dcho. de la lámina

vel=1./sqrt(mu0*epsilon(x1,x2,x,epsvacio,epsilon1));

deltat=cfl*dx./vel;

        epsmat=epsilon(x1,x2,x,epsvacio,epsilon1);
        sigmat=sigma(x1,x2,x,sigma1);
        cEn=dt./(dx*(epsmat+0.5*sigmat.*dt));
        cHn=dt./mu0/dx;
        
dem=length(deltat);        

% Al hacerlo así, la lámina queda más o menos centrada.

%% Graficar la lámina

fun1=@(ds) k1;
fun2=@(ds) k2;
funy=@(ds) ds;
ds=linspace(-1,1,1000);
equis1=fun1(ds);
equis2=fun2(ds);
yies=funy(ds);

%% Performs time integration.
tic;
for t=0:dt:finalTime
    % --- Updates E field ---
    for i=2:cells
%         epsmat=epsilon(x1,x2,x(i),epsvacio,epsilon1);
%         sigmat=sigma(x1,x2,x(i),sigma1);
%         cEn=dt/(dx*(epsmat+sigmat*0.5*dt));
        ez(i,2)=((epsmat(i)-sigmat(i)*0.5*deltat(i))/(epsmat(i)+0.5*sigmat(i)*deltat(i)))*ez(i,1)+cEn(i)*(hy(i-1,1)-hy(i,1));
    end
%     ez(2:end,2) = ez(2:end,1) + cE.* (hy(1:(end-1),1)-hy(2:end,1));
    
    % --- Sources ---
%     ez(excPoint,2) = ez(excPoint,2) + exp(- 0.5*((t-delay)/spread)^2);
    
    % --- Boundary conditions ---
    ez(    1, 2) = -ez(      2, 2);
    ez(cells, 2) = -ez(cells-1, 2); % PEC
    
    % --- Updates H field ---
    for i=1:cells-1
       % cHn=deltat(i)/mu0/dx
        hy(i,2)=hy(i,1)+cHn*(ez(i,2)-ez(i+1,2));
    end
%     hy(1:(end-1),2) = hy(1:(end-1),2) + cH.* (ez(1:(end-1),2) - ez(2:end,2));
    
    ez(:,1)=ez(:,2);
    hy(:,1)=hy(:,2);
    
    % --- Output requests ---
    subplot(2,1,1);
    hold off;
    plot(ez(:,2));
%     hold on
%     plot(equis1,yies,'r',equis2,yies,'r');
%     hold off
    hold on;
    axis([0 cells -1 1]);
    title(sprintf('FDTD Time = %.2f nsec',t*1e9))
    subplot(2,1,2);
    hold off;
    plot(hy(:,2));
%     hold on
%     plot(equis1,yies,'r',equis2,yies,'r');
%     hold off
    hold on;
    axis([0 cells -0.005 0.005]);
    
    drawnow;
end
toc;