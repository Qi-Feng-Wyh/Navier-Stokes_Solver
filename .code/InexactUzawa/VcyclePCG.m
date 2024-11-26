function [XU,XV]=VcyclePCG(F,G,N,tau,eps,error2,v1,v2,maxk)
    k=0;
    XU=zeros(N+1,N);
    XV=zeros(N,N+1);
    PU=zeros(N+1,N);
    PV=zeros(N,N+1);
    RU=F;
    RV=G;
    rho=sum(RU.*RU,'all')+sum(RV.*RV,'all');
    norm_b=sqrt(norm(F,2)^2+norm(G,2)^2);
    norm_r=sqrt(rho);                   
    maxd=log2(N)-1;
    while norm_r/norm_b>eps && norm_r>tau*error2 && k<maxk
        k=k+1;
        [ZU,ZV]=VCycle(RU,RV,N,maxd,v1,v2,zeros(N+1,N),zeros(N,N+1));
        if k==1
            PU=ZU;
            PV=ZV;
            rho=sum(RU.*ZU,'all')+sum(RV.*ZV,'all');
        else
            rho_hat=rho;
            rho=sum(RU.*ZU,'all')+sum(RV.*ZV,'all');
            beta=rho/rho_hat;
            PU=ZU+beta*PU;
            PV=ZV+beta*PV;
        end
        [WU,WV]=apply_A(PU,PV,N);
        alpha=rho/(sum(PU.*WU,'all')+sum(PV.*WV,'all'));
        XU=XU+alpha*PU;
        XV=XV+alpha*PV;
        RU=RU-alpha*WU;
        RV=RV-alpha*WV;
        norm_r=sqrt(norm(RU,2)^2+norm(RV,2)^2);
    end
    
    