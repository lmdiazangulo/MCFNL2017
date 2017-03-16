function [E] = analyticalGaussianMagnetic(x,tOr,L,spread)

    constants;
    
    t = mod(tOr, 2*L/c0);
    
    x0 = L/2;
    
    right1 =  0.5*exp(-(x-x0-c0*t).^2/2/spread^2);
    left1  = -0.5*exp(-(x-x0+c0*t).^2/2/spread^2);
    
    right2 = -0.5*exp(-(x-x0-L+c0*t).^2/2/spread^2);
    left2  = +0.5*exp(-(x-x0+L-c0*t).^2/2/spread^2);
    
    right3 =  0.5*exp(-(x-x0+2*L-c0*t).^2/2/spread^2);
    left3  = -0.5*exp(-(x-x0-2*L+c0*t).^2/2/spread^2);
    
    E = (right1 + left1 + right2 + left2 + right3 + left3)/eta0;
end