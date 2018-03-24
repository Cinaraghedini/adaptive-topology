
function [v] =  compute_eigVector(A,normalized)

D = degree(A);
v=[];

L  =  D - A;

if normalized
    L = ((D)^(-1/2)) * L * ((D)^(-1/2));
end
[eVectors,eValues]=eig(L);
[~,Y]=sort(diag(eValues));
v=eVectors(:,Y(2));