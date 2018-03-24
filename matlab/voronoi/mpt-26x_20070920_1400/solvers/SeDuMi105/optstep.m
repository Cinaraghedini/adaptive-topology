%             [x,y] = optstep(A,b,c, y0,y,d,v,dxmdz, K,L,symLden,...
%                        dense,Ablkjc,Aord,ADA,DAt, feasratio, R,pars)
% OPTSTEP Implements Mehrotra-Ye type optimality projection for
%  IPM-LP solver.
%
% SEE ALSO sedumi
% **********  INTERNAL FUNCTION OF SEDUMI **********
function [x,y] = optstep(A,b,c, y0,y,d,v,dxmdz, K,L,symLden,...
         dense,Ablkjc,Aord,ADA,DAt, feasratio, R,pars)
 %  
 %   This file is part of SeDuMi 1.05
 %   Copyright (C) 2001 Jos F. Sturm
 %     Dept. Econometrics & O.R., Tilburg University, the Netherlands.
 %     Supported by the Netherlands Organization for Scientific Research (NWO).
 %   Affiliation SeDuMi 1.03 and 1.04Beta (2000):
 %     Dept. Quantitative Economics, Maastricht University, the Netherlands.
 %   Affiliations up to SeDuMi 1.02 (AUG1998):
 %     CRL, McMaster University, Canada.
 %     Supported by the Netherlands Organization for Scientific Research (NWO).
 % 
 %   This program is free software; you can redistribute it and/or modify
 %   it under the terms of the GNU General Public License as published by
 %   the Free Software Foundation; either version 2 of the License, or
 %   (at your option) any later version.
 % 
 %   This program is distributed in the hope that it will be useful,
 %   but WITHOUT ANY WARRANTY; without even the implied warranty of
 %   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 %   GNU General Public License for more details.
 % 
 %   You should have received a copy of the GNU General Public License
 %   along with this program; if not, write to the Free Software
 %   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 %

% ------------------------------------------------------------
% For LP-problems, we PROJECT onto OPTIMAL FACE (x0 > 0)  *OR*
%   onto a direction (x0 = 0) with the same sign for c'*x and b'*y.
% ------------------------------------------------------------
 if abs(abs(feasratio)-1) < 0.1
   x0 = sqrt(d.l(1)) * v(1);
   z0 = x0 / d.l(1);
   deptol = 1E-10 * max(x0,z0);
   if (feasratio < -0.5) & (x0 < z0*z0)
     x0 = 0;                          % Try project onto direction.
   end
% ------------------------------------------------------------
% Set d(Non-basic-LP) = 0
% Nonbasic is where dx < dz i.e. dxmdz < 0, ideally: dx_N = -v_N, dz_N = 0.
% ------------------------------------------------------------
   lpNB = find(dxmdz < 0);
   d.l(lpNB) = 0;
   x = sqrt(d.l) .* v;
% --------------------------------------------------
% Compute ADA, with d[N] = 0. Hence A[B]*D[B]^2*A[B]'.
% --------------------------------------------------
   DAt = getDAtm(A,Ablkjc,dense,DAt.denq,d,K);
   ADA = getada1(ADA, A,Ablkjc(:,3),Aord.lqperm, d, K.qblkstart);
   ADA = getada2(ADA, DAt,Aord, K);
   udsqr = invcholfac(d.u,K, d.perm);
   [ADA,absd] = getada3(ADA, A,Ablkjc(:,3),Aord,udsqr,K);
% ------------------------------------------------------------
% Block Sparse Cholesky: ADA(L.perm,L.perm) = L.L*diag(L.d)*L.L'
% ------------------------------------------------------------
   L = sparchol(L,ADA,pars.chol,absd);
% ------------------------------------------------------------
% Factor dense columns
% ------------------------------------------------------------
   [Lden,L.d] = deninfac(symLden, L,dense,DAt,d,absd,K.qblkstart,pars.chol);
% ------------------------------------------------------------
% Solve ADAt*psi = -x0*b+A*D*v, dx = v-D*At*psi.  LEAST SQUARES.
% ------------------------------------------------------------
   [psi,dx,err.kcg,err.b] = wrapPcg(L,Lden,A,dense,d, DAt,K,...
     (-x0) * b,v, pars.cg,pars.eps / pars.cg.restol);
   x = sqrt(d.l) .* dx;
% ----------------------------------------
% CHECK WHETHER x[B] >= 0 AND WHETHER RESIDUAL DID NOT DETERIORATE.
% ----------------------------------------
   if (min(x) < 0.0) | ...
      (norm(err.b,inf) > 5 * max(max(y0,1e-10 * x0) * R.maxb, y0 * R.maxRb))
      x = [];  % Incorrect guess of LP-basis
      return
   else
% ==================== DUAL PART (LP ONLY) ====================
% ------------------------------------------------------------
% Solve A[B]'*dy = x0*c[B]-A[B]'*y, so that zB=0. We solve
% the equivalent system obtained after pre-multiplying with A[B]*P(d[b]).
% ------------------------------------------------------------
     rhs = sqrt(d.l) .* (x0*c - Amul(A,dense,y,1));
% ------------------------------------------------------------
% Solve ADA*dy = A*D*rhs.   THIS IS DEBATABLE !
% ------------------------------------------------------------
     dy = wrapPcg(L,Lden,A,dense,d, DAt,K,...
       zeros(length(b),1),rhs, pars.cg,pars.eps / pars.cg.restol);
     y = y+dy;
% ------------------------------------------------------------
% CHECK WHETHER Z[N] >= 0 AND RESID ON Z[B]=0 DID NOT DETERIORATE.
% ALSO, we need z0 >= 0 and x0+z0 > 0.
% ------------------------------------------------------------
     z = x0*c - Amul(A,dense,y,1);
     z(1) = 0;
     zB = z; zB(lpNB) = 0; normzB = norm(zB,inf);
     cx = c'*x; by = b'*y;
     z0 = by - cx;
     if (min(z(lpNB)) < 0.0) | ...
       normzB > 5 * max(1E-10 * (x0+(x0==0)) * norm(c), y0 * norm(R.c))
       x = [];  % Incorrect guess of LP-basis
       return
     end
     if x0 == 0
       if z0 <= 0
	 x = [];  % z0 <= 0 --> P/D Direction not improving anymore
         return
       end
     elseif z0 < - (5E-8) * (1+abs(by) +max(abs(b)))
       x = [];        % x0>0, z0 << 0 then not optimal.
       return
     end
   end
 else
% ------------------------------------------------------------
% If optimality step impossible , then return emptyset.
% ------------------------------------------------------------
   x = []; y = [];
 end
