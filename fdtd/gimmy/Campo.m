t=0:1e-12:100e-9;
E0=1;
s=0.1;
x0=0.5;
x=0.95;
Ey=analyticalGaussianRight(x,t,E0,x0,s);
plot(t,Ey,'r')