function [mat] = loadMatLib(omega)
Nfreq = numel(omega);
% Loads constants.
constants;
% ==== Creates materials ==================================================
mat = struct('relPermittivity',1 , ...
             'relPermeability',1);
% -------- Defines material 1 ---------------------------------------------
mat(1).dispersive           =   0;
mat(1).electricConductivity = 0.0;
mat(1).relPermittivity      =   1;
mat(1).relPermeability      =   1;
mat(1).Npoles               =   0;

