function [ D ] = DmatA(N,dim)
% function [ D ] = DmatA(n,d)
% Creates D matrix analitically for order n and dimension d
% Calculated from Silvester Ferrari's book p.134
nId = nodeIndices(N,dim);

switch dim
    case 1
        Np = N+1;
    case 2
        Np = (N+1)*(N+2)/2;
    case 3
        Np = (N+1)*(N+2)*(N+3)/6;
end

D=zeros(Np,Np);
for m=1:Np
        switch dim
            case 1
                a = alphaPol1D( Rpol(nId(m,1),N), Rpol(nId(m,2),N) );
                d = zeros(size(a,1)-1,size(a,2)); % d will store the partial derivative of a w.r.t simplex 1.
                for j=1:size(d,1)
                    d(j,:)=j*a(j+1,:);
                end
                for k=1:Np
                    z1=nId(k,1)./N; z2=nId(k,2)./N;
                    D(m,k)=0;
                    for n=1:size(d,1);
                        for p=1:size(d,2);
                            D(m,k)= D(m,k)+ d(n,p)*(z1)^(n-1)*z2^(p-1);
                        end
                    end
                 end     
                   
            case 2
                a = alphaPol2D( Rpol(nId(m,1),N), Rpol(nId(m,2),N), ...
                    Rpol(nId(m,3),N));
                d = zeros(size(a,1)-1,size(a,2),size(a,3));
                for j=1:size(d,1)
                    d(j,:,:)=j*a(j+1,:,:);
                end
                for k=1:Np
                    z1=nId(k,1)./N; z2=nId(k,2)./N; z3=nId(k,3)./N;
                    D(m,k)=0;
                    for n=1:size(d,1);
                        for p=1:size(d,2);
                            for q=1:size(d,3);
                                D(m,k)= D(m,k)+ d(n,p,q) * z1^(n-1) * z2^(p-1) * z3^(q-1);
                            end
                        end
                    end
                 end   
            case 3
                a = alphaPol3D( Rpol(nId(m,1),N), Rpol(nId(m,2),N), ...
                    Rpol(nId(m,3),N), Rpol(nId(m,4),N) );
                d = zeros(size(a,1)-1,size(a,2),size(a,3),size(a,4));
                for j=1:size(d,1)
                    d(j,:,:,:)=j*a(j+1,:,:,:);
                end
                for k=1:Np
                    z1=nId(k,1)./N; z2=nId(k,2)./N; z3=nId(k,3)./N; z4=nId(k,4)./N;
                    D(m,k)=0;
                    for n=1:size(d,1);
                        for p=1:size(d,2);
                            for q=1:size(d,3);
                                for r=1:size(d,4);
                                    D(m,k)= D(m,k)+ d(n,p,q,r)*(z1)^(n-1)*z2^(p-1)*z3^(q-1)*z4^(r-1);
                                end
                            end
                        end
                    end
                 end  
        end 
end