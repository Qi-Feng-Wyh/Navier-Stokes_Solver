
function C=curl_BT(U,V,N)
C=zeros(N,N);
h=1/N;
for i=1:N
    for j=1:N
        C(i,j)=-1/h*(U(i+1,j)-U(i,j)+V(i,j+1)-V(i,j));
    end
end