% [At,b,c,K] = rls(P,q)
% Creates dual standard form for robust least squares problem "Pu=q".
function [At,b,c,K] = rls(P,q)

[m, n] = size(P);
% ---------- minimize y_1 + y_2 ----------
b = -sparse([1; 1; zeros(n,1)]);
% ---------- (y_1, q - P y_3) in Qcone ----------
At = sparse([-1, zeros(1,1+n); ...
             zeros(m,2), P]);
c = [0;q];
K.q = [1+m];
% ---------- (y_2, (1,y_3)) in Qcone ----------
At = [At; 0, -1, zeros(1,n); ...
          zeros(1,2+n); ...
          zeros(n,2), -eye(n)];
c = sparse([c; 0;1;zeros(n,1)]);
K.q = [K.q, 2+n];
