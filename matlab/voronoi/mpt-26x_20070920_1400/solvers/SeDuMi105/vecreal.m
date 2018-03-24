% r = vecreal(xpr,xpi,cpx,K)
% Puts both "real" xpr and "imaginary" xpi into r. No sign changes.
% NB in c'*x = Re(c)'*Re(x) + Im(c)'*Im(x), the sign will indeed be OK.
% **********  INTERNAL FUNCTION OF SEDUMI **********
function r = vecreal(xpr,xpi,cpx,K)

 %  
 %   This file is part of SeDuMi 1.02   (03AUG1998)
 %   Copyright (C) 1998 Jos F. Sturm
 %   CRL, McMaster University, Canada.
 %   Supported by the Netherlands Organization for Scientific Research (NWO).
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

 N = length(xpr) + cpx.dim;
 if issparse(xpr)
   r = sparse([],[],[],N,1,nnz(xpr) + nnz(xpi));
 else
   r = zeros(N,1);
 end
% ----------------------------------------
% LP: only real
% ----------------------------------------
 r(1:K.l) = xpr(1:K.l);
% ----------------------------------------
% LORENTZ: complex (x1, x2) -->  (x1, [real(x2);imag(x2)])
% ----------------------------------------
 firstk = K.l + 1;
 rfirstk = firstk;
 k = 1;
 for knz = 1:length(cpx.q)
   newk = cpx.q(knz);
   len = sum(K.q(k:newk));
   r(rfirstk:rfirstk+len-1) = xpr(firstk:firstk+len-1);  % real part
   rfirstk = rfirstk + len;
   nk = K.q(newk);
   firstk = firstk + len;
   r(rfirstk:rfirstk+nk-2) = xpi(firstk-nk+1:firstk-1);  % imag part
   rfirstk = rfirstk + nk - 1;
   k = newk + 1;
 end
 if k <= length(K.q)
   len = sum(K.q(k:end));
   r(rfirstk:rfirstk+len-1) = xpr(firstk:firstk+len-1);   % remaining blocks
   rfirstk = rfirstk + len;
   firstk = firstk + len;
 end
% ----------------------------------------
% PSD: first real symmetric blocks, i.e. skip Hermitian blocks.
% ----------------------------------------
 hfirstk = firstk;    % save firstk for later use
 k = 1;
 for knz = 1:length(cpx.s)
   newk = cpx.s(knz);
   if newk > k
     len = sum(K.s(k:newk-1).^2);
     r(rfirstk:rfirstk+len-1) = xpr(firstk:firstk+len-1);  % copy real part
     rfirstk = rfirstk + len;
     firstk = firstk + len;
   end
   firstk = firstk + K.s(newk)^2;         % skip over complex block.
   k = newk + 1;
 end
 if k <= length(K.s)
   len = sum(K.s(k:end).^2);
   r(rfirstk:rfirstk+len-1) = xpr(firstk:firstk+len-1);   % remaining blocks
   rfirstk = rfirstk + len;
   firstk = firstk + len;
 end
% ----------------------------------------
% PSD: complex Hermitian part
% ----------------------------------------
 firstk = hfirstk;    % restore to start of PSD variables
 k = 1;
 for knz = 1:length(cpx.s)
   newk = cpx.s(knz);
   if newk > k
     firstk = firstk + sum(K.s(k:newk-1).^2);   % skip over real blocks
   end
   nksqr = K.s(newk)^2;
   r(rfirstk:rfirstk+nksqr-1) = xpr(firstk:firstk+nksqr-1);   % real part
   rfirstk = rfirstk + nksqr;
   r(rfirstk:rfirstk+nksqr-1) = xpi(firstk:firstk+nksqr-1);   % imag part
   rfirstk = rfirstk + nksqr;
   firstk = firstk + nksqr;
   k = newk + 1;
 end
