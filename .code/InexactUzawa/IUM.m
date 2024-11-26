
function [U,V,P]=IUM(F,G,N,alpha,tau,eps,v1,v2,maxk)
    U=zeros(N+1,N);
    V=zeros(N,N+1);
    P=zeros(N,N);    
    norm_r0=sqrt(sum(F.*F,'all')^2+sum(G.*G,'all'));
    norm_rh=norm_r0;
    while norm_rh/norm_r0>=eps
        curl=curl_BT(U,V,N);
        error2=sum(curl.*curl,'all');
        [U1,V1]=apply_B(P,N);
        [U,V]=VcyclePCG(F-U1,G-V1,N,tau,eps,error2,v1,v2,maxk);
        P=P+alpha*curl_BT(U,V,N);
        [AU,AV]=apply_A(U,V,N);
        [BPU,BPV]=apply_B(P,N);
        RU=F-AU-BPU;
        RV=G-AV-BPV;
        D=curl_BT(U,V,N);
        norm_rh=sqrt(sum(RU.*RU,'all')^2+sum(RV.*RV,'all')+sum(D.*D,'all'));
    end
    
    