function [theta_hat, P]=adaptation(lambda,theta_hat_prev, P_prev, phi, f)
gamma=(1/lambda)*P_prev;
P=gamma*(eye(2)-phi*((eye(1)+phi'*gamma*phi)^-1)*phi'*gamma);
e0=f-phi'*theta_hat_prev;
theta_hat=theta_hat_prev+P*phi*e0;
end