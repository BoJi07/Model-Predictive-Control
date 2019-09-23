clc, clear, close all

Me=0.025;Be=0.1;
A=[-Me/Be 0;1 0];
B=[1/Me;0];
C=[0 1];
D=0;

h=0.001;
[Phi,Gamma,Cd,Dd]=c2dm(A,B,C,D,h,'zoh');

K=place(Phi, Gamma, [0.1 0.11]);
eigs(Phi-Gamma*K)

sysd=ss(Phi-Gamma*K,Gamma,Cd,Dd,h);
[y,t]=step(sysd,1);
N=1/y(end);