% [At,b,c,K] = toepest(P)
% Creates dual standard form for Toeplitz-covariance estimation
function [At,b,c,K] = toepest(P)

m = size(P,1);
% ---------- maximize y(m+1) ----------
b = [sparse(m,1); 1];
% ----- Let e be all-1, and allocate space for the A-matrix -----
e = ones(m,1);
K.q = [1 + m*(m+1)/2];
K.xcomplex = 2:K.q(1);      %Norm-bound entries are complex valued
At = sparse([],[],[],K.q(1) + m^2,m+1,1 + 2*m^2);
% ---------- constraints ----------
% -y(m+1) >= norm( vec(P) - sum(y_i * Ti) )     (Qcone)
% sum(y_i * Ti) is psd                          (Scone)
% ---------------------------------
At(:,1) = [sparse(2:(m+1),1,1,K.q(1),1); -vec(speye(m))];
c = [0; diag(P)];
firstk = m+2;
for k = 1:(m-1)
  lastk = firstk + m-k-1;
  Ti = spdiags(e,k,m,m);
  At(:,k+1) = [sqrt(2) * sparse(firstk:lastk,1,1,K.q(1),1); -2*vec(Ti)];
  c = [c; sqrt(2) * diag(P,k)];
  firstk = lastk + 1;
end
At(:,m+1) = [1; sparse(K.q(1) + m^2-1,1)];   % "objective" variable y(m+1)
c = [c; zeros(m^2,1)];              % all-0 in the psd-part
K.s = [m];
K.scomplex=1;                       %Complex Hermitian PSD
% ---------- y(2:m) complex, y(1) and y(m+1) real ----------
K.ycomplex = 2:m;