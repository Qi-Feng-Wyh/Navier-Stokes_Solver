
function [U1,V1,P1]=Restrict(U,V,P,N)
M=N/2;
U1=zeros(M+1,M);
V1=zeros(M,M+1);
P1=zeros(M,M);
for i=1:M
    for j=1:M
        P1(i,j)=1/4*(P(2*i,2*j)+P(2*i-1,2*j)+P(2*i,2*j-1)+P(2*i-1,2*j-1));
    end
end
for j=1:M
    for i=2:M
        U1(i,j)=1/4*(U(2*i-1,2*j)+U(2*i-1,2*j-1))+1/8*(U(2*i-2,2*j)+U(2*i,2*j)+U(2*i-2,2*j-1)+U(2*i,2*j-1));
    end
    U1(1,j)=1/2*(U(1,2*j)+U(1,2*j-1));
    U1(M+1,j)=1/2*(U(N+1,2*j)+U(N+1,2*j-1));
end
for i=1:M
    for j=2:M
        V1(i,j)=1/4*(V(2*i,2*j-1)+V(2*i-1,2*j-1))+1/8*(V(2*i,2*j-2)+V(2*i-1,2*j-2)+V(2*i,2*j)+V(2*i-1,2*j));
    end
    V1(i,1)=1/2*(V(2*i,1)+V(2*i-1,1));
    V1(i,M+1)=1/2*(V(2*i,N+1)+V(2*i-1,N+1));
end