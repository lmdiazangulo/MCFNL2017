function [ x ] = setNodes1D(N,vx)
% function [ x ] = setNodes1D(N,vx)
% Sets N+1 nodes in equispaced positions in the vertices indicated by vx.

K = size(vx,1);

x = zeros(N+1,K);
% Sets nodes in equispaced positions.
for k=1:K
    for i=0:N
        x(i+1,k) = i * (vx(k,2)-vx(k,1))/N + vx(k,1);
    end
end