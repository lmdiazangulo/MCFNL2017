function  [res] = alphaCoeff(N,D)
% function [res] = alphaCoeff(N,D)
% function alphaCoeff(N,D) provides a matrix with as many rows as non-zero
% values has the Lagrange polynomial for order N and dimension D, i,e. the
% number of non zero monomials in each node.
% The number of columns is D+2
% The first column stores the node number that the coefficient is refering
% to.
% Second column stores the constant coefficient for that monomial,
% Following D columns contain the exponent of the simplex coordinates.


nId = nodeIndices(N,D);
switch D
    case 1
        Np = N+1;
    case 2
        Np = (N+1)*(N+2)/2;
    case 3
        Np = (N+1)*(N+2)*(N+3)/6;
end

% Builds the res matrix.
res = [];
for n=1:Np
    
    % Builds polynomials in matrix form,
    switch D
        case 1
            a = alphaPol1D( Rpol(nId(n,1),N), Rpol(nId(n,2),N) );
        case 2
            a = alphaPol2D( Rpol(nId(n,1),N), Rpol(nId(n,2),N), ...
                Rpol(nId(n,3),N));
        case 3
            a = alphaPol3D( Rpol(nId(n,1),N), Rpol(nId(n,2),N), ...
                Rpol(nId(n,3),N), Rpol(nId(n,4),N) );
    end
    
    % Counts the number non zero (nnz) values in that polynomial.
    switch D
        case 1
            for i=1:N+1
                for j=1:N+1
                    if a(i,j) ~= 0
                        res = cat(1,res,[n a(i,j) i-1 j-1]);
                    end
                end
            end
        case 2
            for i=1:N+1
                for j=1:N+1
                    for k=1:N+1
                        if a(i,j,k) ~= 0
                            res = cat(1,res,[n a(i,j,k) i-1 j-1 k-1]);
                        end
                    end
                end
            end
        case 3
            for i=1:N+1
                for j=1:N+1
                    for k=1:N+1
                        for l=1:N+1
                            if a(i,j,k,l) ~= 0
                                res = cat(1,res,[n a(i,j,k,l) i-1 j-1 k-1 l-1]);
                            end
                        end
                    end
                end
            end
    end
            
    
end

