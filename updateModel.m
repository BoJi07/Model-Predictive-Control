function [A,B,C,D]=updateModel(p1,p2)
A=[-p1/p2 0;1 0];
B=[1/p1;0];
C=[0 1];
D=0;
end