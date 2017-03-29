function [ T ] = TmatA(N,D)
% function [ T ] = TmatA(N,D)
% Creates T matrix analitically for orden N and dimension D
nId = nodeIndices(N,D);

switch D
    case 1
        Np = N+1;
    case 2
        Np = (N+1)*(N+2)/2;
    case 3
        Np = (N+1)*(N+2)*(N+3)/6;
end


T = zeros(Np,Np);
for i=1:Np
    for j=1:i      
        switch D
            case 1
                a1 = alphaPol1D( Rpol(nId(i,1),N), Rpol(nId(i,2),N) );
                a2 = alphaPol1D( Rpol(nId(j,1),N), Rpol(nId(j,2),N) );
            case 2
                a1 = alphaPol2D( Rpol(nId(i,1),N), Rpol(nId(i,2),N), ...
                    Rpol(nId(i,3),N));
                a2 = alphaPol2D( Rpol(nId(j,1),N), Rpol(nId(j,2),N), ...
                    Rpol(nId(j,3),N));
            case 3
                a1 = alphaPol3D( Rpol(nId(i,1),N), Rpol(nId(i,2),N), ...
                    Rpol(nId(i,3),N), Rpol(nId(i,4),N) );
                a2 = alphaPol3D( Rpol(nId(j,1),N), Rpol(nId(j,2),N), ...
                    Rpol(nId(j,3),N), Rpol(nId(j,4),N) );
        end
        aMul = convn(a1,a2);
        T(i,j) = integralT(aMul,D);
        T(j,i) = T(i,j);
    end
end