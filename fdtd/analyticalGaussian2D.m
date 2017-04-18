function [E] = analyticalGaussian2D(x,y,x0,y0,spread)

    constants;
    
    E = exp(-( (x-x0).^2 + (y-y0).^2)/2/spread^2);
    

end
