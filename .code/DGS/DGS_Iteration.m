%
function [DGSU,DGSV,DGSP]=DGS_Iteration(U,V,P,F,G,D,N)
U_half=U;
V_half=V;
R=zeros(N,N);
h=1/N;
[BPU,BPV]=apply_B(P,N);
F=F-BPU;
G=G-BPV;
%一步 Gauss-Seidel 迭代
for i=2:N
    for j=1:N
        U_half(i,j)=h^2*F(i,j)+U_half(i-1,j)+U_half(i+1,j);
        if j~=1
            U_half(i,j)=U_half(i,j)+U_half(i,j-1);
        end
        if j~=N
            U_half(i,j)=U_half(i,j)+U_half(i,j+1);
        end
        if j==1||j==N
            U_half(i,j)=1/3*U_half(i,j);
        else
            U_half(i,j)=1/4*U_half(i,j);
        end
    end
end
for i=1:N
    for j=2:N%沟槽的一开始打成N-1了
        V_half(i,j)=h^2*G(i,j)+V_half(i,j-1)+V_half(i,j+1);
        if i~=1
            V_half(i,j)=V_half(i,j)+V_half(i-1,j);
        end
        if i~=N
            V_half(i,j)=V_half(i,j)+V_half(i+1,j);
        end
        if i==1||i==N
            V_half(i,j)=1/3*V_half(i,j);
        else
            V_half(i,j)=1/4*V_half(i,j);
        end
    end
end
%DGS迭代
DGSU=U_half;
DGSV=V_half;
%中间格
for i=2:N-1
    for j=2:N-1
        R(i,j)=-1/h*(DGSU(i+1,j)-DGSU(i,j)+DGSV(i,j+1)-DGSV(i,j))-D(i,j);
        delta=R(i,j)*h/4;
        DGSU(i,j)=DGSU(i,j)-delta;
        DGSU(i+1,j)=DGSU(i+1,j)+delta;
        DGSV(i,j)=DGSV(i,j)-delta;
        DGSV(i,j+1)=DGSV(i,j+1)+delta;
        P(i,j)=P(i,j)+R(i,j);
        P(i+1,j)=P(i+1,j)-1/4*R(i,j);
        P(i-1,j)=P(i-1,j)-1/4*R(i,j);
        P(i,j+1)=P(i,j+1)-1/4*R(i,j);
        P(i,j-1)=P(i,j-1)-1/4*R(i,j);
    end
end
%边缘格
for i=2:N-1
    j=1;
    R(i,j)=-1/h*(DGSU(i+1,j)-DGSU(i,j)+DGSV(i,j+1)-DGSV(i,j))-D(i,j);
    delta=R(i,j)*h/3;
    DGSU(i,j)=DGSU(i,j)-delta;
    DGSU(i+1,j)=DGSU(i+1,j)+delta;
    DGSV(i,j+1)=DGSV(i,j+1)+delta;
    P(i,j)=P(i,j)+R(i,j);
    P(i+1,j)=P(i+1,j)-1/3*R(i,j);
    P(i-1,j)=P(i-1,j)-1/3*R(i,j);
    P(i,j+1)=P(i,j+1)-1/3*R(i,j);
    
    j=N;
    R(i,j)=-1/h*(DGSU(i+1,j)-DGSU(i,j)+DGSV(i,j+1)-DGSV(i,j))-D(i,j);
    delta=R(i,j)*h/3;
    DGSU(i,j)=DGSU(i,j)-delta;
    DGSU(i+1,j)=DGSU(i+1,j)+delta;
    DGSV(i,j)=DGSV(i,j)-delta;
    P(i,j)=P(i,j)+R(i,j);
    P(i+1,j)=P(i+1,j)-1/3*R(i,j);
    P(i-1,j)=P(i-1,j)-1/3*R(i,j);
    P(i,j-1)=P(i,j-1)-1/3*R(i,j);
end

for j=2:N-1
    i=1;
    R(i,j)=-1/h*(DGSU(i+1,j)-DGSU(i,j)+DGSV(i,j+1)-DGSV(i,j))-D(i,j);
    delta=R(i,j)*h/3;
    DGSU(i+1,j)=DGSU(i+1,j)+delta;
    DGSV(i,j)=DGSV(i,j)-delta;
    DGSV(i,j+1)=DGSV(i,j+1)+delta;
    P(i,j)=P(i,j)+R(i,j);
    P(i+1,j)=P(i+1,j)-1/3*R(i,j);
    P(i,j+1)=P(i,j+1)-1/3*R(i,j);
    P(i,j-1)=P(i,j-1)-1/3*R(i,j);
    
    i=N;
    R(i,j)=-1/h*(DGSU(i+1,j)-DGSU(i,j)+DGSV(i,j+1)-DGSV(i,j))-D(i,j);
    delta=R(i,j)*h/3;
    DGSU(i,j)=DGSU(i,j)-delta;
    DGSV(i,j)=DGSV(i,j)-delta;
    DGSV(i,j+1)=DGSV(i,j+1)+delta;
    P(i,j)=P(i,j)+R(i,j);
    P(i-1,j)=P(i-1,j)-1/3*R(i,j);
    P(i,j+1)=P(i,j+1)-1/3*R(i,j);
    P(i,j-1)=P(i,j-1)-1/3*R(i,j);
end
%角落格子
i=1;j=1;
R(i,j)=-1/h*(DGSU(i+1,j)-DGSU(i,j)+DGSV(i,j+1)-DGSV(i,j))-D(i,j);
delta=R(i,j)*h/2;
DGSU(i+1,j)=DGSU(i+1,j)+delta;
DGSV(i,j+1)=DGSV(i,j+1)+delta;
P(i,j)=P(i,j)+R(i,j);
P(i+1,j)=P(i+1,j)-1/2*R(i,j);
P(i,j+1)=P(i,j+1)-1/2*R(i,j);

i=1;j=N;
R(i,j)=-1/h*(DGSU(i+1,j)-DGSU(i,j)+DGSV(i,j+1)-DGSV(i,j))-D(i,j);
delta=R(i,j)*h/2;
DGSU(i+1,j)=DGSU(i+1,j)+delta;
DGSV(i,j)=DGSV(i,j)-delta;
P(i,j)=P(i,j)+R(i,j);
P(i+1,j)=P(i+1,j)-1/2*R(i,j);
P(i,j-1)=P(i,j-1)-1/2*R(i,j);

i=N;j=1;
R(i,j)=-1/h*(DGSU(i+1,j)-DGSU(i,j)+DGSV(i,j+1)-DGSV(i,j))-D(i,j);
delta=R(i,j)*h/2;
DGSU(i,j)=DGSU(i,j)-delta;
DGSV(i,j+1)=DGSV(i,j+1)+delta;
P(i,j)=P(i,j)+R(i,j);
P(i-1,j)=P(i-1,j)-1/2*R(i,j);
P(i,j+1)=P(i,j+1)-1/2*R(i,j);

i=N;j=N;

R(i,j)=-1/h*(DGSU(i+1,j)-DGSU(i,j)+DGSV(i,j+1)-DGSV(i,j))-D(i,j);
delta=R(i,j)*h/2;
DGSU(i,j)=DGSU(i,j)-delta;
DGSV(i,j)=DGSV(i,j)-delta;
P(i,j)=P(i,j)+R(i,j);
P(i-1,j)=P(i-1,j)-1/2*R(i,j);
P(i,j-1)=P(i,j-1)-1/2*R(i,j);
%返回值
DGSP=P;
