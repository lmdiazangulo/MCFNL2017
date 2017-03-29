function [P] = P3Da(N,f)
% This function provides the Permutation matrix which leaves the f vertex
% were originally was the one vertex.

Np = (N+1)*(N+2)*(N+3)/6;


% Defines Q1 as the clockwise rotation taking vertex 1 as axis.
Q1 = zeros(Np,Np);
orNum = 1:Np;
fiNum = zeros(Np,1);
nodesSet = 0;
for i=1:(N+1)
    Nsp = i*(i+1)/2; % Nodes in the slide  
    fiNum((nodesSet+1):(nodesSet+Nsp)) = ...
        Pmatrix2D(i-1)*orNum((nodesSet+1):(nodesSet+Nsp))';   
    nodesSet = nodesSet+Nsp;
end

for i=1:Np
    Q1(orNum(i),fiNum(i))=1;
end

% Defines Q2, same as Q1 but around vertex 2.
orNum=zeros(Np,1);
lastSet=N*(N+1)*(N+2)/6+1;
lastNodeSet=0;
for i=0:N
    for j=1:(i+1)
        orNum(lastNodeSet+j)=lastSet;
        lastSet = lastSet + 1;
        lastNodeSet = lastNodeSet+j;
    end
end

lastSet = N*(N+1)*(N+2)/6+1;
for j=1:N
    temp = orNum;
    for i=(Np):(-1):2
        if temp(i)~=0
            if temp(i-1)==0
                lastSet=lastSet-1;
                orNum(i-1)=lastSet;
            end
        end
    end
end

nodesSet=0; fiNum=zeros(Np,1);
for i=1:(N+1)
    Nsp = i*(i+1)/2; % Nodes in the slide
    temp = Pmatrix2D(i-1)*orNum((nodesSet+1):(nodesSet+Nsp));
    fiNum((nodesSet+1):(nodesSet+Nsp)) = temp;
    nodesSet = nodesSet+Nsp;
end

Q2 = zeros(Np,Np);            
for i=1:Np
    Q2(orNum(i),fiNum(i))=1;
end

% Builds Q3 and Q4.
Q3 = Q2*Q2*Q1*Q1;
Q4 = Q1*Q1*Q2*Q2;


P = zeros(Np,Np);
switch f
    case 1
        P = eye(Np,Np);
    case 2
        P = Q3*Q2*Q2;
    case 3
        P = Q4*Q3*Q3;
    case 4
        P = Q2*Q4*Q4;
end



end