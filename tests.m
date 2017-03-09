clear all;
close all;

constants;

L = 1;
spread = 0.05;

x=(0:0.01:1)';
t=(0:1e-11:10e-9);
E = zeros(size(t,2), size(x,1));
for i=1:length(t)
    E(i,:) = analyticalGaussianRight(x,t(i),L,spread);
end

surf(E);
shading interp;