function [res] = analyticalGaussianRight(x,t,L,spread)

    constants;
    
    x0 = L/2;
    
    right = exp(-(x-x0-c0*t).^2/2/spread^2);
    left  = - exp(-(x-x0-L+c0*t).^2/2/spread^2);
    
    res = right + left;
end