function [ alpha ] = alphaPol1D(r1, r2)

alpha = zeros(numel(r1),numel(r2));
for i=1:numel(r2)
    alpha(:,i) = conv(r1,r2(i));
end

end