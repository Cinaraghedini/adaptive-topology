clear;

disp('2nd problem');
s=1;
c=[-1,-1]';
a=[-2,5;2,-2];
b=[5;1];
ctype=['U','U']';
lb=[0;0]; ub=[];
vartype=['I';'I'];
param.msglev=1;
[xmin,fmin,status,extra]=glpkmex(s,c,a,b,ctype,lb,ub,vartype,param)
pause;

disp('3rd problem');
s=1;
c=[0 0 0 -1 -1]';
a=[-2 0 0 1 0;...
    0 1 0 0 2;...
    0 0 1 3 2];
b=[4 12 18]';
ctype=['S','S','S']';
lb=[0,0,0,0,0]'; ub=[];
vartype=['C','C','C','C','C']';
[xmin,fmin,status,extra]=glpkmex(s,c,a,b,ctype,lb,ub,vartype)
