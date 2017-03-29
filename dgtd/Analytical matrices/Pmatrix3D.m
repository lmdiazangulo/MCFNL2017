function [ P ] = Pmatrix3D(N, sDir)
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
switch sDir
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

% % % Np = (N+1)*(N+2)*(N+3)/6;
% % % P2 = zeros(Np,Np); P3 = zeros(Np,Np);
% % % 
% % % switch N
% % %     case 1
% % %         P2 = [0 1 0 0;
% % %               1 0 0 0;
% % %               0 0 0 1;
% % %               0 0 1 0];
% % %         P3 =  [0 0 1 0;
% % %                0 0 0 1;
% % %                1 0 0 0;
% % %                0 1 0 0];
% % %     case 2
% % %         P2(1,5)=1; P2(2,2)=1;  P2(3,7)=1;
% % %         P2(4,6)=1; P2(8,10)=1; P2(9,9)=1;    
% % %         
% % %         P3(1,8)=1; P3(2,9)=1;  P3(3,3)=1;
% % %         P3(4,6)=1; P3(5,10)=1; P3(7,7)=1; 
% % %           
% % %     case 3
% % %         P2(1,11)=1; P2(2,5)=1; P2(3,13)=1; P2(4,12)=1; P2(6,7)=1;
% % %         P2(8,16)=1;  P2(9,15)=1; P2(10,14)=1; P2(17,20)=1; P2(18,19)=1;
% % %         
% % %         P3(1,17)=1; P3(2,18)=1; P3(3,8)=1; P3(4,14)=1; P3(5,19)=1;
% % %         P3(6,9)=1; P3(7,15)=1; P3(10,12)=1; P3(11,20)=1; P3(13,16)=1;
% % % 
% % %     case 4
% % %         P2(1,21)=1; P2(2,11)=1; P2(3,23)=1; P2(4,22)=1; P2(5,5)=1;
% % %         P2(6,13)=1; P2(7,12)=1; P2(8,26)=1; P2(9,25)=1; P2(10,24)=1;
% % %         P2(14,16)=1; P2(15,15)=1; P2(17,30)=1; P2(18,29)=1; P2(19,28)=1;
% % %         P2(20,27)=1; P2(31,35)=1; P2(32,34)=1; P2(33,33)=1;
% % %         
% % %         P3(1,31)=1; P3(2,32)=1; P3(3,17)=1; P3(4,27)=1; P3(5,33)=1;
% % %         P3(6,18)=1; P3(7,28)=1; P3(8,8)=1; P3(9,14)=1; P3(10,24)=1;
% % %         P3(11,34)=1; P3(12,19)=1; P3(13,29)=1; P3(15,15)=1; P3(16,25)=1;
% % %         P3(20,22)=1; P3(21,35)=1; P3(23,30)=1; P3(26,26)=1;
% % % end
% % % 
% % % for i=1:Np
% % %     for j=1:Np
% % %         P2(j,i) = P2(i,j);
% % %         P3(j,i) = P3(i,j);
% % %     end
% % % end
% % % 
% % % switch sDir
% % %     case 1, P = eye(Np,Np);
% % %     case 2, P = P2;
% % %     case 3, P = P3;
% % %     case 4, P = P2*P3;
% % % end