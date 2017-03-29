function [ e ] = setZeroField(e)
% function [ e ] = setZeroField(e)

K = numel(e);
Np = size(e(1).x, 1);

for k = 1:K
    e(k).E = zeros(Np,1);
    e(k).H = zeros(Np,1);
end
