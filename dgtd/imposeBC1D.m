function [ e ] = imposeBC1D(e, leftCond, rightCond)
% function imposeBC1D(e, leftCond, rightCond)
% Imposes boundary conditions modifying jumps across interfaces.
% e is an array of element's structs.
% right and leftCond must be one of these numbers:
% 0 -> SM-ABC (Silver-Muller Absorbing Boundary Condition).
% 1 -> PEC (Perfect Electric Conductor).
% 2 -> PMC (Perfect Magnetic Codunctor).
constants;
% Preliminar definitions.
K = numel(e);
Np = size(e(1).x, 1);
% Imposes BC over the left boundary.
switch leftCond
    case 0 % SMA-BC
        e(1).dE(1) = e(1).E(1);
        e(1).dH(1) = e(1).H(1);
    case 1 % PEC
        e(1).dE(1) = +2 * e(1).E(1);
        e(1).dH(1) = 0;
    case 2 % PMC
        e(1).dE(1) = 0;
        e(1).dH(1) = +2 * e(1).H(1);
    otherwise
        disp('ERROR: Invalid BC on left boundary.');
end
% Imposes BC over right boundary
switch rightCond
    case 0 % SMA
        % Impedance terms.
%         e(K).rhsD = e(K).impPoles.a .* e(K).D ...
%                     + e(K).impPoles.r .* e(K).H(Np);
%         e(K).rhsB = e(K).admPoles.a .* e(K).B ...
%                     + e(K).admPoles.r .* e(K).E(Np);
        e(K).dE(2) = e(K).E(Np);% 2 * e(K).E(Np) - e(K).Z .* e(K).H(Np) - sum(e(K).rhsD);
        e(K).dH(2) = e(K).H(Np);%e(K).Z .* e(K).H(Np) + sum(e(K).rhsD);
    case 1 % PEC
        e(K).dE(2) = +2 * e(K).E(Np);
        e(K).dH(2) = 0;
    case 2 % PMC
        e(K).dE(2) = 0;
        e(K).dH(2) = +2 * e(K).H(Np);
    otherwise
        disp('ERROR: Invalid BC on right boundary.');
end
