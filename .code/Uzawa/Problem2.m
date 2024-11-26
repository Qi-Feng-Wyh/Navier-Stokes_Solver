N=256;
%创建矩阵，通过键值对的方式
I=zeros(10*N*N-18*N+4,1);
J=zeros(10*N*N-18*N+4,1);
S=zeros(10*N*N-18*N+4,1);
cnt=1;
for i=1:N
    for j=1:N
        if i~=N
            ind_i=j+(i-1)*N;
            if j==1||j==N
                I(cnt)=ind_i;
                J(cnt)=ind_i;
                S(cnt)=3;
                cnt=cnt+1;
            else
                I(cnt)=ind_i;
                J(cnt)=ind_i;
                S(cnt)=4;
                cnt=cnt+1;
            end
            if j>1 %(i,j-1) ok
                ind_j=j-1+(i-1)*N;
                I(cnt)=ind_i;
                J(cnt)=ind_j;
                S(cnt)=-1;
                cnt=cnt+1;
            end
            if j<N %(i,j+1)ok
                ind_j=j+1+(i-1)*N;
                I(cnt)=ind_i;
                J(cnt)=ind_j;
                S(cnt)=-1;
                cnt=cnt+1;
            end
            if i>1 %(i-1,j)ok
                ind_j=j+(i-2)*N;
                I(cnt)=ind_i;
                J(cnt)=ind_j;
                S(cnt)=-1;
                cnt=cnt+1;
            end
            if i<N-1 %(i+1,j)ok
                ind_j=j+i*N;
                I(cnt)=ind_i;
                J(cnt)=ind_j;
                S(cnt)=-1;
                cnt=cnt+1;
            end
        end
        if j~=N
            ind_i=j+(i-1)*(N-1)+N*(N-1);
            if i==1||i==N
                I(cnt)=ind_i;
                J(cnt)=ind_i;
                S(cnt)=3;
                cnt=cnt+1;
            else
                I(cnt)=ind_i;
                J(cnt)=ind_i;
                S(cnt)=4;
                cnt=cnt+1;
            end
            if j>1 %(i,j-1) ok
                ind_j=j-1+(i-1)*(N-1)+N*(N-1);
                I(cnt)=ind_i;
                J(cnt)=ind_j;
                S(cnt)=-1;
                cnt=cnt+1;
            end
            if j<N-1 %(i,j+1)ok
                ind_j=j+1+(i-1)*(N-1)+N*(N-1);
                I(cnt)=ind_i;
                J(cnt)=ind_j;
                S(cnt)=-1;
                cnt=cnt+1;
            end
            if i>1 %(i-1,j)ok
                ind_j=j+(i-2)*(N-1)+N*(N-1);
                I(cnt)=ind_i;
                J(cnt)=ind_j;
                S(cnt)=-1;
                cnt=cnt+1;
            end
            if i<N %(i+1,j)ok
                ind_j=j+i*(N-1)+N*(N-1);
                I(cnt)=ind_i;
                J(cnt)=ind_j;
                S(cnt)=-1;
                cnt=cnt+1;
            end
        end
    end
end
A=sparse(I,J,S);
I=zeros(4*N*(N-1),1);
J=zeros(4*N*(N-1),1);
S=zeros(4*N*(N-1),1);
S=S/N;
cnt=1;
for i=1:N
    for j=1:N
        if i~=N %上半部分，N-1个I_N
            ind_i=j+(i-1)*N;
            ind_j=ind_i+N;
            I(cnt)=ind_i;
            J(cnt)=ind_i;
            S(cnt)=-1;
            cnt=cnt+1;
            I(cnt)=ind_i;
            J(cnt)=ind_j;
            S(cnt)=1;
            cnt=cnt+1;
        end
        if j~=N %下半部分，N 个 Q
            ind_i=j+(i-1)*(N-1)+N*(N-1);
            ind_j=j+(i-1)*N;
            I(cnt)=ind_i;
            J(cnt)=ind_j;
            S(cnt)=-1;
            cnt=cnt+1;
            I(cnt)=ind_i;
            J(cnt)=ind_j+1;
            S(cnt)=1;
            cnt=cnt+1;
        end
    end
end
B=sparse(I,J,S);
C=B';
U=zeros(2*N*(N-1),1);
P=zeros(N*N,1);
F=zeros(2*N*(N-1),1);
%设置初值
for i=1:N
   for j=1:N
       if i~=N
            x=i/N;y=(j-0.5)/N;
            F(j+(i-1)*N)=-4*pi*pi*(2*cos(2*pi*x)-1)*sin(2*pi*y)+x*x;
            if j==1
                F(j+(i-1)*N)=F(j+(i-1)*N)-2*pi*(1-cos(2*pi*x))*N;
            elseif j==N
                F(j+(i-1)*N)=F(j+(i-1)*N)+2*pi*(1-cos(2*pi*x))*N;
            end
       end
       if j~=N
            x=(i-0.5)/N;y=j/N;
            F(N*(N-1)+j+(i-1)*(N-1))=4*pi*pi*(2*cos(2*pi*y)-1)*sin(2*pi*x);
            if i==1
                F(N*(N-1)+j+(i-1)*(N-1))=F(N*(N-1)+j+(i-1)*(N-1))+2*pi*(1-cos(2*pi*y))*N;
            elseif i==N
                F(N*(N-1)+j+(i-1)*(N-1))=F(N*(N-1)+j+(i-1)*(N-1))-2*pi*(1-cos(2*pi*y))*N;
            end
       end
   end
end
F=F/(N*N);
%系数小，收敛快
U=sparse(U);
P=sparse(P);
F=sparse(F);
tic;
while 1
    U=A\(F-B*P);
    P=P+C*U;
    rh=sqrt(norm(A*U+B*P-F,2)*norm(A*U+B*P-F,2)+norm(C*U,2)*norm(C*U,2));
    if rh/norm(F,2)<1e-8
        break;
    end
end
disp(toc);
error=0;
for i=1:N
   for j=1:N
       if i~=N
            x=i/N;y=(j-0.5)/N;
           	tmp=U(j+(i-1)*N)-(1-cos(2*pi*x))*sin(2*pi*y);
            error=error+tmp*tmp;
       end
       if j~=N
            x=(i-0.5)/N;y=j/N;
            tmp=U(N*(N-1)+j+(i-1)*(N-1))+(1-cos(2*pi*y))*sin(2*pi*x);
            error=error+tmp*tmp;
       end
   end
end
error=sqrt(error)/N;
disp(error);