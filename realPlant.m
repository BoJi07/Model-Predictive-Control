function [x0,x]=realPlant(x0_prev,ti,te)
tint=ti;tfin=te;
[tt,state]=ode45(@sys_ode,[tint,tfin],x0_prev);
x0=state(end,:);
x=state(end,:)';
end