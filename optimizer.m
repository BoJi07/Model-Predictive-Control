%% Iterative Optimizer
function U=optimizer(umin,umax,x)
global Z Q W R m rp N

alpha_k=0.005;
beta_k=0.005;
Acon=kron(eye(N),[-eye(m);eye(m)]);
Dg=Acon;
Bcon=kron(ones(N,1),[-umin;umax]);
U=zeros(N*m,1);
mu=zeros(N*m*2,1);
% KKT_test1=0;
% KKT_test2=0;
% KKT_test3=0;
end_test=0;
cnt=0;
while ~end_test  && cnt<=10000
    g=Acon*U-Bcon;
    DJ=-Z'*Q*(rp-W*x-Z*U)+R*U;
    dU=-alpha_k*(DJ+Dg'*mu);
    U=U+dU;
    mu=max(mu+beta_k*g,0);
    %     KKT_test1=isequal(mu,zeros(Np*m*2,1));
    %     KKT_test2=isequal(DJ+Dg'*mu,zeros(Np*m,1));
    %     KKT_test3=mu'*g==0;
    end_test=max(abs(dU))<=1e-3;
    cnt=cnt+1;
end
end