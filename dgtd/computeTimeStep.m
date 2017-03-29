function [ dt ] = computeTimeStep(e, CFL)
% function [ dt ] = computeTimeStep(e, CFL)
% dt must use the minimum distance to be calculated. That is not yet
% implemented and the first element spatial size is used instead.

constants;  
K  = numel(e);      % Number of elements
Np = numel(e(1).x); % Number of nodes per element.

% Computes minimum distance between nodes.
xMin = e(1).x(2) - e(1).x(1);
for k=1:K
    for i=1:(Np-1)
        if ( e(k).x(i+1) - e(k).x(i) ) < xMin
            xMin = e(k).x(i+1) - e(k).x(i) ;
        end
    end
end

dt = CFL*xMin/c0;
