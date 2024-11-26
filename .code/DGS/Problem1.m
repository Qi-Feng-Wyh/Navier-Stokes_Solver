error_N=zeros(1,6);
N=64;
for ind=1:6
    h=1/N;
    [F,G]=init(N);
    D=zeros(N,N);
    [U0,V0,P0]=true_solution(N);
    U=zeros(N+1,N);
    V=zeros(N,N+1);
    P=zeros(N,N);                         
    maxd=log2(N)-1;
    [U,V,P,time,itercnt]=VCycleDGS(U,V,P,F,G,D,N,maxd,5,5);
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
    error_N(ind)=error;
    fprintf('N=%d, error=%e\n', N, error);
    maxd=log2(N)-1;
    [U,V,P,time,itercnt]=VCycleDGS(U,V,P,F,G,D,N,maxd,3,3);
    fprintf('N=%d, bottom=2*2, v1=v2=3, time=%f, itercnt=%d\n', N, time, itercnt);
    [U,V,P,time,itercnt]=VCycleDGS(U,V,P,F,G,D,N,maxd-1,3,3);
    fprintf('N=%d, bottom=4*4, v1=v2=3, time=%f, itercnt=%d\n', N, time, itercnt);
    [U,V,P,time,itercnt]=VCycleDGS(U,V,P,F,G,D,N,maxd,5,5);
    fprintf('N=%d, bottom=2*2, v1=v2=5, time=%f, itercnt=%d\n', N, time, itercnt);
    [U,V,P,time,itercnt]=VCycleDGS(U,V,P,F,G,D,N,maxd-1,5,5);
    fprintf('N=%d, bottom=4*4, v1=v2=5, time=%f, itercnt=%d\n', N, time, itercnt);
    [U,V,P,time,itercnt]=VCycleDGS(U,V,P,F,G,D,N,maxd,7,7);
    fprintf('N=%d, bottom=2*2, v1=v2=7, time=%f, itercnt=%d\n', N, time, itercnt);
    [U,V,P,time,itercnt]=VCycleDGS(U,V,P,F,G,D,N,maxd-1,7,7);
    fprintf('N=%d, bottom=4*4, v1=v2=7, time=%f, itercnt=%d\n', N, time, itercnt);
    N=N*2;
end
x=[6,7,8,9,10,11];
y=log2(error_N);
disp(size(x));
disp(size(y));
scatter(x,y,'markerfacecolor',[0,0,1]);
xlabel('log_2 N');
ylabel('log_2 e_N');
p=polyfit(x,y,1);
fprintf('k=%e\n', p(1));

