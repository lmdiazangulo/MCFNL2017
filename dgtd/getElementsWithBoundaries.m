function [ eL eR ] = getElementsWithBoundaries(e, xBC)

K = numel(e);
Np = numel(e(1).x);
found = 0;
for i=1:K;
    for j = 1:Np
        if (e(i).x(j) == xBC && found == 0)
            eL = i;
            if (i == K)
                eR = 0;
            else
                eR = i + 1;
            end
            found = 1;
        end
    end
end
