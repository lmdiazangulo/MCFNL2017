%%%%
% Laura Martin, 18th of April
%
% With this program one can simulate electromagnetic waves in 
% one dimension. The algorithm is based on fdtd(2,2)[red line] 
% and fdtd(2,4)[blue line].
% 
% To use fdtd(2,4) a cfl<=6/7 is necassary. 
% 
% Two scenarios: 
% [Un/-comment the respective lines in the code.]
% (1.) A wave hits an infinitesimal small PEC layer. In case of 4th-order a small 
% part of the wave will be transmitted. This is from the physical point of view 
% not posible.
% (2.) Comparing 2nd and 4th order simulation with the analytical result. 
%
%
% Reference list: 
% -- Wilson, Joshua, et al. 
%    "An accurate and stable fourth order finite difference time domain method." 
%    Microwave Symposium Digest, 2008 IEEE MTT-S International. IEEE, 2008.
% -- Georgakopoulos, Stavros V., et al. 
%    "Higher-order finite-difference schemes for electromagnetic radiation, 
%     scattering, and penetration. 1. theory." 
%    IEEE Antennas and Propagation Magazine 44.1 (2002): 134-142.
% -- Yefet, Amir, and Peter G. Petropoulos. 
%    "A staggered fourth-order accurate explicit finite difference scheme 
%     for the time-domain Maxwell's equations." 
%    Journal of Computational Physics 168.2 (2001): 286-315.
%%%%

close all;
clear variables;

%% LOAD CONSTANTS
constants;

%% SPACIAL PARAMETERS
L =10;                 
dx = 0.025;
x = (0:dx:L)'; 
cells = size(x,1)+2;  

%% TEMPORAL PARAMETERS
cfl = 4/7;            % proportion of time and space resolution (2nd cfl=1)
dt = cfl*dx/c0;       % resolution of time
T = L/c0*0.8;         % temporal limit

%% PARAMS
cE2 = dt/(eps0*dx);
cE4 = dt/(eps0*dx*24);
cH2 = dt/(mu0*dx);
cH4 = dt/(mu0*dx*24);

%% STORAGE FOR SIMULATION
ez2 = zeros(cells  ,2); % 2nd-order resolution
hy2 = zeros(cells-1,2);
ez4 = zeros(cells  ,2);
hy4 = zeros(cells-1,2); % 4th-order resolution

%% (1.) INITIAL E-FIELD  (run with lines 120,121 and 138->end)
%spread = 0.5;
%ez2(2:cells-1,1) = analyticalGaussian(x,-dt/2,L,spread); 
%ez4(2:cells-1,1) = analyticalGaussian(x,-dt/2,L,spread); 

%% (2.) SOURCE PARAM   (run with lines 76,77,88,89 and 120)
loc = floor(length(x)/4);
spread = 1e-9;
delay = 5e-9;
scaPoint = floor(length(x)*3/4);

%% SIMULATION LOOP
%tic; 
for t=0:dt:(T+dt/2)
  %% UPDATE E-FIELD
  for i=3:(cells-2)
    % 2th order %% hy/ez
    ez2(i,2) = ez2(i,1) + cE2*(hy2(i,1)-hy2(i-1,1));
    % 4th order %% hy/ez
    ez4(i,2) = ez4(i,1) + ...
               cE4* (hy4(i-2,1)-27*hy4(i-1,1)+27*hy4(i,1)-hy4(i+1,1));
  end
  %% HANDLE E-FIELD SOURCE 
  %% (1.)
  ez2(loc,2) = ez2(loc,2) + exp(-0.5*((t-delay)/spread)^2);
  ez4(loc,2) = ez4(loc,2) + exp(-0.5*((t-delay)/spread)^2);
  %% E-FIELD BOUNDARIES
  %-- PEC 2th Order --%
  ez2(      2,2) = 0;  
  ez2(cells-1,2) = 0; 
  %-- PEC 4th Order -- Dirichlet Condition --%
  ez4(     1 ,2) = 0; 
  ez4(     2 ,2) = 0;        
  ez4(cells-1,2) = 0; 
  ez4(  cells,2) = 0; 
  %-- infinitesimal PEC layer in the middle --% 
  %% (1.)
  ez2(   201,2 ) = 0; 
  ez4(   201,2 ) = 0; 
  
  %% UPDATE H-FIELD
  for i=2:(cells-2)
    % 2th order %% hy/ez
    hy2(i,2) = hy2(i,1) + cH2*(ez2(i+1,2)-ez2(i,2));
    % 4th order %% hy/ez  
    hy4(i,2) = hy4(i,1) + ...
               cH4* (ez4(i-1,2)-27*ez4(i,2)+27*ez4(i+1,2)-ez4(i+2,2));
  end
  %% H-FIELD BOUNDARIES
  %-- PEC 4th order -- Neumann Condition --%
  hy4(      1,2) = hy4(2,2);
  hy4(cells-1,2) = hy4(cells-2,2);
  %% HANDLE H-FIELD SOURCE
  
  %% PREPARE FOR NEXT TIMESTEP
  ez2(:,1) = ez2(:,2);
  hy2(:,1) = hy2(:,2);
  ez4(:,1) = ez4(:,2);
  hy4(:,1) = hy4(:,2);
  
  
  %% PLOT
  %if isequal(mod(t,40*dt),0)
    subplot(2,1,1);
    hold off;
    plot(x, ez2(2:end-1,1), '-r' ); 
    hold on;
    plot(x, ez4(2:end-1,1), '--b' ); 
    hold on;
    %% (1.)
    line([x(200) x(200)],[-1 1]);
    %% (2.)
%    xAn = 0:0.001:L;
%    plot(xAn,analyticalGaussian(xAn,t+dt/2,L,spread),'-g');
    axis( [x(1) x(end) -1 1] );
    title(sprintf('FDTD Step = %d of %d',floor(t/dt),floor(T/dt)))  

    subplot(2,1,2);
    hold off;
    plot(x(1:end-1), hy2(2:end-1,1), '-r' );
    hold on;
    plot(x(1:end-1), hy4(2:end-1,1), '--b' );
    axis( [x(1) x(end) -0.005 0.005] );
    
    drawnow; % show plot
  %end
end
%toc;

%% SHOW ANALYTICAL RESULT AND ERROR
%% (2.)
%Eanalytical = analyticalGaussian(x,t+dt/2,L,spread);
%fprintf('For dx = %e \n', dx);
%fprintf('   cfl = %.3f \n', cfl);
%fprintf('     T = %e \n', T);
%fprintf('2nd order:  L^2 error: %e\n', ...
%     sum(abs(ez2(2:cells-1,1)-Eanalytical))/cells);
%fprintf('4th order:  L^2 error: %e\n', ...
%     sum(abs(ez4(2:cells-1,1)-Eanalytical))/cells);
%%subplot(1,1,1);
%%plot(x,ez2(:,1), '-r' );
%%hold on;
%%plot(x, ez4(:,1), '..b' );
%%hold on;
%%plot(x(2:cells-1), Eanalytical, '.g' );
%%drawnow;