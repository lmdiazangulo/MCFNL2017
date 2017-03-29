function [ alpha ] = alphaPol3D(r1, r2, r3, r4)

alpha = zeros(numel(r1),numel(r2),numel(r3),numel(r4));
for i=1:numel(r2)
    for j=1:numel(r3)
        for k=1:numel(r4)
            alpha(:,i,j,k) =conv( conv( conv(r1,r2(i)), r3(j)), r4(k) );
        end
    end
end

end
