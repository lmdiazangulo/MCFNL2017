function [ e ] = computeJumps(e)
% function computeJumps(e)

Np = size(e(1).E,1);
K = numel(e);

% k=1 is a particular case.
e(1).dE(2) = e(1).E(Np) - e(2).E(1);
e(1).dH(2) = e(1).H(Np) - e(2).H(1);
% Bulk of elements.
for k = 2:(K-1)
    e(k).dE(1) = e(k).E(1)  - e(k-1).E(Np);
    e(k).dH(1) = e(k).H(1)  - e(k-1).H(Np);
    e(k).dE(2) = e(k).E(Np) - e(k+1).E(1);
    e(k).dH(2) = e(k).H(Np) - e(k+1).H(1);
end
% k=K is another particular case.
e(K).dE(1) = e(K).E(1) - e(K-1).E(Np);
e(K).dH(1) = e(K).H(1) - e(K-1).H(Np);
