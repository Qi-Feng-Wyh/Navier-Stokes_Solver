%
function [AU,AV]=apply_A(U,V,N)
AU=zeros(N+1,N);
AV=zeros(N,N+1);
h=1/N;
%迭代水平速度
for i=2:N
    for j=2:N-1
        AU(i,j)=1/(h^2)*(4*U(i,j)-U(i,j-1)-U(i-1,j)-U(i+1,j)-U(i,j+1));
    end
    AU(i,1)=1/(h^2)*(3*U(i,1)-U(i-1,1)-U(i+1,1)-U(i,2));
    AU(i,N)=1/(h^2)*(3*U(i,N)-U(i-1,N)-U(i+1,N)-U(i,N-1));
end
for j=1:N
    AU(1,j)=U(1,j);
    AU(N+1,j)=U(N+1,j);
end
%迭代竖直速度
for j=2:N
    for i=2:N-1
        AV(i,j)=1/(h^2)*(4*V(i,j)-V(i,j-1)-V(i-1,j)-V(i,j+1)-V(i+1,j));
    end
    AV(1,j)=1/(h^2)*(3*V(1,j)-V(1,j-1)-V(1,j+1)-V(2,j));
    AV(N,j)=1/(h^2)*(3*V(N,j)-V(N,j-1)-V(N,j+1)-V(N-1,j));
end
for i=1:N
    AV(i,1)=V(i,1);
    AV(i,N+1)=V(i,N+1);
end