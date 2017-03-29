% +=======================================================================+
% | 1-Dimensional Discontinuous Galerkin for Maxwell Equations.           |
% | Created by Luis Manuel Diaz Angulo for the University of Granada.     |
% +-----------------------------------------------------------------------+
% | PURPOSE:                                                              |
% | This script simulates propagation of an EM wave in 1-D using a DGTD   |
% | algorithm.                                                            |
% +=======================================================================+
clear all;
% ------- Physical and mathematical constants -----------------------------
constants;
%% ============== Initialization ==========================================
% -------Simulation general parameters -----------------------------------
N = 5;           % N    := Degree of polynomials forming the basis.
h = 0.1;     % h    := Spatial resolution. [meters]
xMin = 0;        % xMin := Left boundary position.  [meters]
xMax = 1;      % xMax := Right boundary position. [meters]
finalTime = 8e-9; % Final time [seconds]
CFL = 0.1;       % CFL condition for time integration stability. 
framesShown = 5; % Plots fields every "framesShown" time steps.
alpha = 0.0;     % 1 = Central , 0 = Upwind
z =    0.0;      % Probe position in the dispersive media.
% -------------------------------------------------------------------------
vx = [xMin:h:(xMax-h); (xMin+h):h:xMax]';   % vx := Element vertices.
K = numel(vx)/2;      % K    := Number of finite elements.
Np = N+1;             % Np   := Number of nodes per element.
x = setNodes1D(N,vx); % x    := Nodes in equispaced positions.
% Creates element structure.
for k=1:K
    e(k) = struct('x', x(:,k));
end
% ------- Material properties ---------------------------------------------
Ns = 1e3;
freqMin = 1e2;
freqMax = 1e9;
fq = linspace(freqMin, freqMax, Ns);
omega = 2 * pi * fq;
matLib = loadMatLib(omega');
m1 = matLib(1);
e = setMaterial(e, m1, 1, K);
% ------- Initial conditions ----------------------------------------------
[ e ] = setZeroField(e); % Sets initial fields to zero.
% ------- Analytical matrices definitions ---------------------------------
addpath 'Analytical matrices';
T    = TmatA(N,1);     % General Mass Matrix for that order.
D{1} = Dmatrix1D(N,1); % Derivative matrix for simplex coordinate 1.
D{2} = Dmatrix1D(N,2); % Derivative matrix for simplex coordiante 2.
% ------- Computes element parameters -------------------------------------
for k=1:K
    M = (e(k).x(Np)-e(k).x(1)) /2 .* T; % M := Mass Matrix.
    e(k).invM = inv( M );                % invM := Inverse Mass Matrix
    e(k).S =  (D{2} - D{1})' / (e(k).x(Np)-e(k).x(1));  % S := Stiffness Matrix.
    e(k).LIFT = e(k).invM ;              % LIFT := Fluxes operator.
end
% ------------- Time parameters initialization ----------------------------
dt     = computeTimeStep(e, CFL); % dt     := Time step.
Nsteps = floor(finalTime / dt);   % Nsteps := Number of iterations needed. 
time   = 0:dt:(Nsteps*dt);        % time   := Time vector.
% ------------- Initialization of aux. variables --------------------------
for k=1:K
    % Right hand side variables.
    e(k).rhsE = zeros(Np,1);
    e(k).rhsH = zeros(Np,1);
    % Jumps between elements.
    e(k).dE = zeros(2,1);
    e(k).dH = zeros(2,1);
    % Fluxes.
    e(k).fluxE = zeros(Np,1);
    e(k).fluxH = zeros(Np,1);
    % Residual fields.
    e(k).resE = zeros(Np,1);
    e(k).resH = zeros(Np,1);
end
% Incident field.
As = 1;
Ag = 1;
freq = 0;
spread = 3e-10;
disp   = 20e-10;
[eS, eT] = getElementsWithBoundaries(e, 0.1);
[eT2, eS2] = getElementsWithBoundaries(e, 0.9);
%% ============= Time integration ========================================= 
% ------------- Main iterative loop ---------------------------------------
figure(3);
for tstep = 1:Nsteps
    % ------------ Computes five stage Runge Kutta time integration -------
    for INTRK = 1:5
        % ---------------- RHS terms computation --------------------------
        e = computeJumps(e);
        e = imposeBC1D(e, SMABC, SMABC);
        timeDelayed = time(tstep) - e(eS).x(Np) / c0;
        e = setIncidentField(e, eT, eS, ... 
             incidentFieldAmplitude(As, freq, timeDelayed, ...
                                    Ag, spread, disp));
	    timeDelayed = time(tstep) - e(eT2).x(Np) / c0;
        e = setIncidentField(e, eT2, eS2, ... 
             incidentFieldAmplitude(As, freq, timeDelayed, ...
                                    Ag, spread, disp));
        e = computeFluxes(e, alpha);
        e = computeRHS(e);
        for k = 1:K
            % ---------------- Updates fields -----------------------------
            e(k).resE = rk4a(INTRK) * e(k).resE + dt * e(k).rhsE;
            e(k).resH = rk4a(INTRK) * e(k).resH + dt * e(k).rhsH;
            e(k).E    = e(k).E + rk4b(INTRK) * e(k).resE;
            e(k).H    = e(k).H + rk4b(INTRK) * e(k).resH;
        end
    end
    % ------------ Plots fields -------------------------------------------
    if mod(tstep, framesShown) == 0
        plotFields(e, time(tstep));
    end
end
