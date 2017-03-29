function [ e ] = computeRHS(e)

K = numel(e);

constants;

for k=1:K
    e(k).rhsE = ( - e(k).S * e(k).H + e(k).LIFT * e(k).fluxE) ./ e(k).elecPermittivity;
    e(k).rhsH = ( - e(k).S * e(k).E + e(k).LIFT * e(k).fluxH) ./ e(k).magnPermeability;
end

