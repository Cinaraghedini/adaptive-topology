clear;

disp('1st problem');
s=-1;
c=[10,6,4]';
a=[1,1,1;...
   10,4,5;...
   2,2,6];
b=[100,600,300]';
ctype=['U','U','U']';
lb=[0,0,0]';
ub=[];
vartype=['C','C','C']';
param.msglev=3;
lpsolver=1;
save_pb=0;
[xmin,fmin,status,extra]=glpkmex(s,c,a,b,ctype,lb,ub,vartype,param,lpsolver,save_pb)
