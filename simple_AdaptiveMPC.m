clc, clear, close all

global Me_act Be_act uk
global Z Q W R m rp N
Me_lim=[0.025 0.085];
Be_lim=[0.1 0.35];
u_lim=[-5 5];
sampling_freq=500;
Q_=10000;
R_=0.0001;
dt=1/sampling_freq;
Tfinal=4;
N=20;
Time=0:dt:5;
ref=0.1*sin(Time*pi);
x0=[0;0];

Me_est=0.055;Be_est=0.225;
Me_act=0.055;Be_act=0.225;
Me_act_=linspace(Me_lim(1),Me_lim(2),2000);
Be_act_=linspace(Be_lim(1),Be_lim(2),2000);

[A,B,C,D]=updateModel(Me_est,Be_est);
[n,m,p,W,Z,Q,R,Kr]=getMPCMatrices(A,B,C,D,N,dt,Q_,R_);

x(:,1)=x0;
y(:,1)=C*x(:,1);
t(1)=0;
e(:,1)=y(:,1)-ref(1);
theta_hat=[1/Me_est;Be_est/Me_est];
lambda=0.995;
P=100*eye(2);
save_theta_hat(:,1)=theta_hat;
save_theta(:,1)=[1/Me_act;Be_act/Me_act];
for k=1:Tfinal/dt
    rp=(ref(k+1:k+N))';
    %u(:,k)=Kr*(R+Z'*Q*Z)^(-1)*Z'*Q*(rp-W*x(:,k));
    u(:,k)=Kr*optimizer(u_lim(1),u_lim(2),x(:,k));
    if u(:,k)>u_lim(2)
        u(:,k)=u_lim(2);
    end
    if u(:,k)<u_lim(1)
        u(:,k)=u_lim(1);
    end
    uk=u(:,k);
    
    [x0, x(:,k+1)]=realPlant(x0,t(k),t(k)+dt);
    
    ydt=x(1,k)+normrnd(0,0.1);
    yddt=1/Me_act*uk-Be_act/Me_act*x(1,k)+normrnd(0,0.1);
    
%     [theta_hat, P]=adaptation(lambda, theta_hat, P, [uk -ydt]', yddt);
%     [A,B,C,D]=updateModel(Me_est,Be_est);
%     [n,m,p,W,Z,Q,R,Kr]=getMPCMatrices(A,B,C,D,N,dt,Q_,R_);
    
%     Me_act=Me_act_(k);
%     Be_act=Be_act_(k);
%     Me_act=(Me_lim(2)-Me_lim(1)).*rand()+ Me_lim(1);
%     Be_act=(Be_lim(2)-Be_lim(1)).*rand()+ Be_lim(1);
    save_theta_hat(:,k+1)=theta_hat;
    save_theta(:,k+1)=[1/Me_act;Be_act/Me_act];

    j(k)=1/2*Q_*e(:,k)^2+1/2*R_*u(:,k)^2;
    y(:,k+1)=C*x(:,k+1);
    e(:,k+1)=y(:,k+1)-ref(k+1);
    u(:,k+1)=0;
    t(k+1)=t(k)+dt;
    clc
    k
end


j(k+1)=1/2*Q_*e(:,k+1)^2+1/2*R_*u(:,k+1)^2;
sum(e.*e)

figure(1)
subplot(2,2,1)
plot(t,y,t,ref(1:Tfinal/dt+1))
title('Trajectory Following')
legend('output','reference')
xlabel('Time(s)'),ylabel('lateral displacement')

subplot(2,2,2)
plot(t,e)
title('Tracking Error')
xlabel('Time(s)'),ylabel('e')

subplot(2,2,3)
stairs(t,u)
title('Control Input')
xlabel('Time(s)'),ylabel('u')

subplot(2,2,4)
plot(t,j)
title('Instantaneous Cost')
xlabel('Time(s)'),ylabel('J')
str=['\Sigma = ',num2str(sum(j))];
text(1,max(j),str);

figure(2)
plot(t,save_theta_hat,t,save_theta)

% figure(2)
% subplot(1,2,1)
% plot(t,theta(1,:))
% xlabel('Time(s)'),ylabel('Me')
% 
% subplot(1,2,2)
% plot(t,theta(2,:))
% xlabel('Time(s)'),ylabel('Be')