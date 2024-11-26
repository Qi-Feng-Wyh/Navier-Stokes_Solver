
function [U_true,V_true,P_true]=true_solution(N)
U_true=zeros(N+1,N);
V_true=zeros(N,N+1);
P_true=zeros(N,N);
h=1/N;
for i=1:N+1
    for j=1:N
        x=(i-1)*h;
        y=(j-0.5)*h;
        U_true(i,j)=(1-cos(2*pi*x))*sin(2*pi*y);
    end
end
for j=1:N+1
    for i=1:N
        x=(i-0.5)*h;
        y=(j-1)*h;
        V_true(i,j)=-(1-cos(2*pi*y))*sin(2*pi*x);
    end
end
for i=1:N
    for j=1:N
        x=(i-0.5)*h;
        P_true(i,j)=x^3/3-1/12;
    end
end