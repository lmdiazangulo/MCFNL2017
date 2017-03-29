function [ Tnode ] = integralT(aMul, D)
% Performs the integral for a pair of nodes.
% aMul is the hyper-matrix of coefficients of the integrand: x1 rows, x2 columns,
% ...
% D is dimension

Tnode = 0;

switch D
    case 1
        for i=1:size(aMul,1)
            for j=1:size(aMul,2)
                if aMul(i,j)~= 0
                    coeffIntegral = factorial(i-1)*factorial(j-1) ...
                        /factorial(i+j+1-2);
                    Tnode = Tnode + aMul(i,j)*coeffIntegral;
                end
            end
        end

    case 2
        for i=1:size(aMul,1)
            for j=1:size(aMul,2)
                for k=1:size(aMul,3)
                    if aMul(i,j,k)~= 0
                        coeffIntegral= factorial(i-1)*factorial(j-1)* ...
                            factorial(k-1)*factorial(2) / ...
                            factorial(i+j+k+2-3);
                        Tnode = Tnode+aMul(i,j,k)*coeffIntegral;
                    end
                end
            end
        end

    case 3
        for i=1:size(aMul,1)
            for j=1:size(aMul,2)
                for k=1:size(aMul,3)
                    for l=1:size(aMul,4)
                        if aMul(i,j,k,l)~=0
                            coeffIntegral= factorial(i-1)*factorial(j-1)* ...
                                factorial(k-1)*factorial(l-1)*factorial(3) / ...
                                factorial(i+j+k+l+3-4);
                            Tnode = Tnode+aMul(i,j,k,l)*coeffIntegral;
                        end
                    end
                end
            end
        end

end