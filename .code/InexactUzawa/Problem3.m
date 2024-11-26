alpha=0.5;
tau=0.01;
eps=1e-8;
N=512;
maxd=log2(N)-1;
h=1/N;
[F,G]=init(N);
[U0,V0,P0]=true_solution(N);
[U,V,P]=IUM(F,G,N,alpha,tau,eps,100,100,2*N*(N-1));
error=0;
for i=2:N
    for j=1:N
        error=error+(U(i,j)-U0(i,j))^2;
    end
end
for j=2:N
    for i=1:N
        error=error+(V(i,j)-V0(i,j))^2;
    end
end
error=h*sqrt(error);
disp(error);