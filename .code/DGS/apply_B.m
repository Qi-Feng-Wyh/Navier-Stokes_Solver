%
function [BPU,BPV]=apply_B(P,N)
BPU=zeros(N+1,N);
BPV=zeros(N,N+1);
h=1/N;
for i=2:N
    for j=1:N
        BPU(i,j)=1/h*(P(i,j)-P(i-1,j));
    end
end
for j=2:N
    for i=1:N
        BPV(i,j)=1/h*(P(i,j)-P(i,j-1));
    end
end