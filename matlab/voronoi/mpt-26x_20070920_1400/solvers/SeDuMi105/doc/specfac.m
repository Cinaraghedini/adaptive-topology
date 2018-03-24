% [At,b,c,K] = specfac(b)
% Creates primal standard form for minimal phase spectral factorization.
function [At,b,c,K] = specfac(b)

m = length(b);
% ---------- minimize sum (m-i)*X(i,i) ----------
c = vec(spdiags((m-1:-1:0)',0,m,m));
% ----- Let e be all-1, and allocate space for the A-matrix -----
e = ones(m,1);
At = sparse([],[],[],m^2,m,m*(m+1)/2);
% ---------- sum(diag(X,k)) = b(k) ----------
for k = 1:m
  At(:,k) = vec(spdiags(e,k-1,m,m));
end
K.s = [m];
