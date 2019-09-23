function xdot = sys_ode(t,x)
global Me_act Be_act uk
xdot=[(-Be_act/Me_act)*x(1)+(1/Me_act)*uk;
            x(1);];
end
