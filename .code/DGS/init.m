%
function [F,G]=init(N)
F=zeros(N+1,N);
G=zeros(N,N+1);
h=1/N;
for i=2:N
    for j=1:N
        x=(i-1)*h;
        y=(j-0.5)*h;
        F(i,j)=-4*pi^2*(2*cos(2*pi*x)-1)*sin(2*pi*y)+x^2;
    end
    F(i,1)=F(i,1)+1/h*(2*pi*(cos(2*pi*x)-1));
    F(i,N)=F(i,N)+1/h*(-2*pi*(cos(2*pi*x)-1));
end
for j=2:N
    for i=1:N
        x=(i-0.5)*h;
        y=(j-1)*h;
        G(i,j)=4*pi^2*(2*cos(2*pi*y)-1)*sin(2*pi*x);
    end
    G(1,j)=G(1,j)+1/h*(-2*pi*(cos(2*pi*y)-1));
    G(N,j)=G(N,j)+1/h*(2*pi*(cos(2*pi*y)-1));
end
   