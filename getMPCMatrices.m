function [n,m,p,W,Z,Q,R,Kr]=getMPCMatrices(A,B,C,D,Np,dt,Q_,R_)
h=dt;
[Phi,Gamma,Cd,Dd]=c2dm(A,B,C,D,h,'zoh');
n=size(Phi,1);m=size(Gamma,2);p=size(Cd,1);

W=[];
Z=[];
Kr=[];
for np=1:Np
    W=cat(1,W,C*Phi^np);
    Zrow=[];
    for npp=1:Np
        if npp<=np
            Zrow=cat(2,Zrow,C*(Phi^(np-npp))*Gamma);
        else
            Zrow=cat(2,Zrow,zeros(p,m));
        end
    end
    Z=cat(1,Z,Zrow);
    if np==1
        Kr=cat(2,Kr,eye(m));
    else
        Kr=cat(2,Kr,zeros(m));
    end
end
Q=Q_*eye(Np*p);
R=R_*eye(Np*m);
end