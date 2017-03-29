function [ e ] = setMaterial(e, m, firstElement, lastElement)
% =========================================================================
% function [ e ] = setMaterial(e, m, firstElement, lastElement)
% e  = element structure to be modified.
% m  = material to be assigned. This material has to be created using the 
%      function createMatLib.
% first/lastElement = first and last element where material must be
%                     included.
% =========================================================================
constants;
for k = firstElement:lastElement
    e(k).elecPermittivity = m.relPermittivity * eps0;
    e(k).magnPermeability = m.relPermeability * mu0;
    e(k).Z = sqrt( e(k).magnPermeability ./ ...
                   e(k).elecPermittivity );
    e(k).Y = 1 / e(k).Z;
end
