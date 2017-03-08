clear;
clc;
%% Dibujo de analytical gaussian right

% Definimos los parámetros

s=0.1;
x0=0.5;
E0=1;
x=0.95;
t=0:1e-12:100e-9;

% El campo es

Campo=analyticalGaussianRight(E0,x0,s,x,t);

% Hacemos el plot

plot(t,Campo)
xlabel('t')
ylabel('Ey')