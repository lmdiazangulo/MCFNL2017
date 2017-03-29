function [P]=Pmatrix1D(N,sDir)
Np = N+1;

switch sDir
    case 1
        P=eye(Np);
    case 2
        P=zeros(Np,Np);
        for i=1:Np
            P(i,Np-i+1)=1;
        end
end