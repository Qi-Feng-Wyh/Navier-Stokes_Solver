
function [SOLU,SOLV,SOLP,time,itercnt]=VCycleDGS(U0,V0,P0,F0,G0,D0,N0,maxd,v1,v2)
U=cell(1,maxd+3);
V=cell(1,maxd+3);
P=cell(1,maxd+3);
F=cell(1,maxd+3);
G=cell(1,maxd+3);
D=cell(1,maxd+3);
RF=cell(1,maxd+3);
RG=cell(1,maxd+3);
RD=cell(1,maxd+3);
r0=sqrt(norm(F0,'fro')*norm(F0,'fro')+norm(G0,'fro')*norm(G0,'fro'));
d=1;
U{d}=U0;
V{d}=V0;
P{d}=P0;
F{d}=F0;
G{d}=G0;
D{d}=D0;
N=N0;
tic;
for i=1:v1
    [U{d},V{d},P{d}]=DGS_Iteration(U{d},V{d},P{d},F{d},G{d},D{d},N);
end
[AU,AV]=apply_A(U{d},V{d},N);
[PU,PV]=apply_B(P{d},N);
RF{d}=F{d}-AU-PU;
RG{d}=G{d}-AV-PV;
RD{d}=D{d}-curl_BT(U{d},V{d},N);%一开始减号两边写反了
itercnt=0;
while 1
    itercnt=itercnt+1;
    while d<maxd
        [F{d+1},G{d+1},D{d+1}]=Restrict(RF{d},RG{d},RD{d},N);
        N=N/2;
        d=d+1;
        U{d}=zeros(N+1,N);
        V{d}=zeros(N,N+1);
        P{d}=zeros(N,N);
        for i=1:v1
            [U{d},V{d},P{d}]=DGS_Iteration(U{d},V{d},P{d},F{d},G{d},D{d},N);
        end
        [AU,AV]=apply_A(U{d},V{d},N);
        [PU,PV]=apply_B(P{d},N);
        RF{d}=F{d}-AU-PU;
        RG{d}=G{d}-AV-PV;
        RD{d}=D{d}-curl_BT(U{d},V{d},N);
    end
    while d>1
        [U1,V1,P1]=Improve(U{d},V{d},P{d},N);
        d=d-1;
        N=N*2;
        U{d}=U{d}+U1;
        V{d}=V{d}+V1;
        P{d}=P{d}+P1;
        for i=1:v2
            [U{d},V{d},P{d}]=DGS_Iteration(U{d},V{d},P{d},F{d},G{d},D{d},N);
        end
    end
    [AU,AV]=apply_A(U{d},V{d},N);
    [PU,PV]=apply_B(P{d},N);
    RF{d}=F{d}-AU-PU;
    RG{d}=G{d}-AV-PV;
    RD{d}=D{d}-curl_BT(U{d},V{d},N);
    rh=sqrt(norm(RF{d},'fro')*norm(RF{d},'fro')+norm(RG{d},'fro')*norm(RG{d},'fro')+norm(RD{d},'fro')*norm(RD{d},'fro'));
    if(rh/r0<=1e-8)
        break;
    end
end
time=toc;
SOLU=U{1};
SOLV=V{1};
SOLP=P{1};
