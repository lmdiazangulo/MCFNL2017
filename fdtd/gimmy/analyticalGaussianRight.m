function [Ey]=analyticalGaussianRight(x,t,E0,L,s)
constants;

%Ey0=E0*exp(-(x-x0).^2/((2*s));
%H0=Ey0/c0;
%Ey=E0*exp(-(c0*t-(x-L/2)/(2*s)).^2);
Ey=E0*exp(-(-c0*t+x-L/2/(2*s)).^2)-E0*exp(-(c0*t+x-3*L/2/(2*s)).^2);

end
