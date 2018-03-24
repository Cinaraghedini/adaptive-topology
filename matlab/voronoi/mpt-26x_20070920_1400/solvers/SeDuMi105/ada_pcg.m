% [y,r,k, DAy] = ada_pcg(L,Lden,At,dense,d, DAt,K, b, cgpars, y0)
%
% Solve y from AP(d)A' * y = b
% using PCG-method and Cholesky L as conditioner.
% If L is sufficiently accurate, then only 1 CG-step is needed.
% In general, proceeds until ||DAy - DA'(AD^2A')^{-1}b||
% has converged (to zero). k = #(CG-iterations).
%
function [y,r,k, DAy] = ada_pcg(L,Lden,At,dense,d, DAt,K, b, cgpars, y0)
% --------------------------------------------------
% INITIALIZE:
%  y=0 with r = b.
% --------------------------------------------------
 restol = y0 * cgpars.restol;
 r = b;
 normrmin = norm(r,inf);
 k = 0;
 STOP = 0;
 finew = 0;
 y = [];
 ymin = [];
% --------------------------------------------------
% Conjugate Gradient until convergence
% --------------------------------------------------
 while ~STOP
% --------------------------------------------------
% P r e - C o n d i t i o n i n g:
% Solve L*Lr = r, L*THETA*tmp = r.
% First iterate: p solves L*THETA*L'*p = r.
% --------------------------------------------------
   Lr = fwdpr1(Lden,sparfwslv(L, r));
   tmp = Lr ./ L.d;
   if (k == 0) | cgpars.rescor
     ssqrNew = Lr'*tmp;
     p = sparbwslv(L, bwdpr1(Lden,tmp));
   else
% --------------------------------------------------
% Iterate k > 0: make p conjugate to previous iterate(s):
% --------------------------------------------------
     ssqrOld = ssqrNew;
     ssqrNew = Lr'*tmp;
     p = (ssqrNew/ssqrOld) * p;
     p = p + sparbwslv(L, bwdpr1(Lden,tmp));
   end
% --------------------------------------------------
% SCALING OPERATION AND MATRIX*VECTOR.
% Let DDAp = P(d)*A'*p and ssqrDAp = ||P(d)^{1/2}*A'*p||^2.
% Set alpha = ssqrNew / ssqrDAp
% --------------------------------------------------
   Ap = vecsym(Amul(At,dense,p,1), K);
   [DDAp, DApq, DAps, ssqrDAp] = PopK(d,Ap,K);
   if ssqrDAp > 0.0
     k = k + 1;
% --------------------------------------------------
% Take step:  y := y + alpha*p
%--------------------------------------------------
     if cgpars.rescor
       alpha = 1;
     else
       alpha = ssqrNew / ssqrDAp;
     end
     if ~isempty(y)
       if isstruct(y)
	 [y.hi,y.lo] = quadadd(y.hi,y.lo,alpha*p);
       else
	 y = y + alpha * p;
       end
     elseif cgpars.qprec > 0
       y.hi = alpha * p; y.lo = zeros(length(p),1);
     else
       y = alpha * p;     
     end
% --------------------------------------------------
% Update residual r := r - alpha * A*[P(d)*Ap]. MATRIX*VECTOR.
% --------------------------------------------------
     tmp = Amul(At,dense,DDAp) + DAt.q'*DApq(:,1) +...
         DAt.denq*DApq(dense.qs,1);
     r = r - alpha * tmp;
% --------------------------------------------------
% Convergence check (HEURISTIC)
% --------------------------------------------------
     fiprev = finew;
     finew = (b+r)'*y;
     normr = norm(r,inf);
     if normr < normrmin
       ymin = y;
       normrmin = normr;
     end
     if normr < restol
       STOP = 1;
     elseif (finew-fiprev < cgpars.stagtol * fiprev) & ~cgpars.rescor
       STOP = 2;
     elseif k >= cgpars.maxiter
       STOP = 2;
     end
   else
     fprintf('Warning: DAp = 0 in PCG\n');
     STOP = 1;         % If DAp == 0 then can't go on.
   end
 end
% --------------------------------------------------
% OUTPUT: return projection D * At*y on request.
% --------------------------------------------------
 y = ymin;
 if isempty(y)
    y = zeros(length(b),1);    % done nothing
 end
 if nargout >= 2
   if k == 1
     DAy = alpha*[sqrt(d.l).*Ap(1:K.l); asmDxq(d,Ap,K,DApq); DAps];
   else
     if isstruct(y)
       Ap = vecsym(Amul(At,dense,y.hi,1), K);
     else
       Ap = vecsym(Amul(At,dense,y,1), K);
     end
     [DDAp, DApq, DAps, ssqrDAp] = PopK(d,Ap,K);
     r = b-Amul(At,dense,DDAp) - DAt.q'*DApq(:,1)...
	 - DAt.denq*DApq(dense.qs,1);
     DAy = [sqrt(d.l).*Ap(1:K.l); asmDxq(d,Ap,K,DApq); DAps];
     if isstruct(y)
       Ap = vecsym(Amul(At,dense,y.lo,1), K);
       [DDAp, DApq, DAps, ssqrDAp] = PopK(d,Ap,K);
       r = r-Amul(At,dense,DDAp) - DAt.q'*DApq(:,1)...
	   - DAt.denq*DApq(dense.qs,1);
       DAy = DAy + [sqrt(d.l).*Ap(1:K.l); asmDxq(d,Ap,K,DApq); DAps];
     end
   end
 end
