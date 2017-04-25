close all;
clear variables;

constants; % Loads constants

%% Problem definition.
x = (0:0.05:10)';
finalTime = 200e-9;

% Materials.

% Boundary conditions.
tic
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
S=zeros(size(x,1),3);
if (exist('initialEz','var'))
    ez(:,1) = initialEz(:);
end
ezLinear=ez;
hyLinear=hy;

%Initial Frequency 
InitialFreq=abs(fftn(ez(:,1)));

% Inits Absorbing Boundary Conditions.
exm2=0;             
exm1=0;
exn2=0;
exn1=0;

% No Linear constants
omegaR=150e6;
alpha=0.5;
tau=sqrt(2)/omegaR;
deltaR=1/tau;
% strength of the third-order nonlinearity
X0=3e-2;
NLC=X0*alpha;

% Determines recursion coefficients
cfl = 0.9;
dx = sum(x(2:end)-x(1:(end-1)))/(length(x)-1);
dt = cfl*dx/c0;

% No Linear constants used in time iteration 
c1=(2-(omegaR*dt)^2)/(deltaR*dt+1);
c2=(deltaR*dt-1)/(deltaR*dt+1);
c3=((1-alpha)*X0*(omegaR*dt)^2)/(deltaR*dt+1); %no linear coef

cE = dt/eps0/dx;
cH = dt/mu0/dx;

%% Performs time integration.
%tic;
for t=0:dt:finalTime
    % --- Updates E field ---
    for i=2:cells
        % Linear part
        ezLinear(i,2)=ezLinear(i,1)+cE*(hyLinear(i-1,1)-hyLinear(i,1));
        
        % NonLinear part
        S(i,3)=c1*S(i,2)+c2*S(i,1)+c3*ez(i,1);
        c=cE*(hy(i-1,1)-hy(i,1))+ez(i,1)*(1+S(i,2))+NLC*ez(i,1)^3;
       %fun=@(y) NLC*y.^3+y-c;
       %x0=ez(i,1);
       %ez(i,2)=fsolve(fun,x0,options);
       C=[NLC 0 1+S(i,3) -c];
       a=roots(C);       
       ez(i,2)=a(3);
    end
%     ez(2:end,2) = ez(2:end,1) + cE.* (hy(1:(end-1),1)-hy(2:end,1));
    
    % --- Sources ---
%     ez(excPoint,2) = ez(excPoint,2) + exp(- 0.5*((t-delay)/spread)^2);
    
    % --- Boundary conditions ---
    ez(    1, 2) = 0;
    ez(cells, 2) = 0; % PEC
    
    % --- Boundary conditions ---Linear 
    ezLinear(    1, 2) = 0;
    ezLinear(cells, 2) = 0; % PEC
    
    
     
     
    % --- Updates H field ---
    for i=1:cells-1
        hy(i,2)=hy(i,1)+cH*(ez(i,2)-ez(i+1,2));
        hyLinear(i,2)=hyLinear(i,1)+cH*(ezLinear(i,2)-ezLinear(i+1,2));
    end
%     hy(1:(end-1),2) = hy(1:(end-1),2) + cH.* (ez(1:(end-1),2) - ez(2:end,2));

    ezLinear(:,1)=ezLinear(:,2);
    hyLinear(:,1)=hyLinear(:,2);
    
    ez(:,1)=ez(:,2);
    hy(:,1)=hy(:,2);
    S(:,1)=S(:,2);
    S(:,2)=S(:,3);
    
    
    %toc
%     figure(1);
%     subplot(2,1,1);
%     hold off;
%     plot(ez(:,2));
%     hold on;
%     plot(ezLinear(:,2),'r');
%     axis([0 cells -1 1]);
%     title(sprintf('FDTD Time = %.2f nsec',t*1e9))
%     subplot(2,1,2);
%     hold off;
%     plot(hy(:,2));
%     hold on;
%     plot(hyLinear(:,2),'r');
%     axis([0 cells -0.01 0.01]);
    

end
% --- Output requests ---
%Initial Frequency 
figure(2)
hold on
FinalFreqLinear= abs(fftn(ezLinear(:,2)));
FinalFreqNonLinear=abs((fftn(ez(:,2))));
plot(InitialFreq,'b');
plot(abs(FinalFreqLinear),'r');
plot(abs(FinalFreqNonLinear),'g');
h = legend('InitialFreq','FinalLinear','FinalNonLinear',3);
set(h,'Interpreter','none')
hold off
    
toc;

initialEz = exp(- (x-5).^2);
syms x;
f = exp(-(x-5)^2);
F=fourier(f)
