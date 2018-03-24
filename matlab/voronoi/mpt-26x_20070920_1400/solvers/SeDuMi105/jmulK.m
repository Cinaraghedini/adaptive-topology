% z = jmulK(x,y, K)  for full x,y, or
% z = jmulK(x,y, K,f) with x in spectral space and Lorentz frame f.
% z = x jmul y (Jordan multiplication).
% **********  INTERNAL FUNCTION OF SEDUMI **********

function z = jmulK(x,y, K,f)

% THE M-FILE VERSION OF THIS FUNCTION IS HERE ONLY AS ILLUSTRATION.
% SEE THE C-SOURCE FOR THE MEX-VERSION.

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

error('At OS prompt, type "make" to create SeDuMi mex-files.')

 z = zeros(length(y),1);
% ----------------------------------------
% LP: z = x.*y
% ----------------------------------------
 z(1:K.l) = x(1:K.l).*y(1:K.l);
 firstk = K.l + 1;
 nextx = K.l+1; nextf = 1;
% ----------------------------------------
% LORENTZ:
% z = x jmul y
% ----------------------------------------
 if(length(x) < length(y))
   for k=1:length(K.q)
     nk = K.q(k); lastk = firstk + nk - 1;
     z(firstk:lastk) = qjmul(x(nextx:nextx+1), y(firstk:lastk), ...
       f(nextf:nextf+nk-2));
     firstk = lastk + 1; nextf = nextf + nk-1; nextx = nextx + 2;
   end
 else
   for k=1:length(K.q)
     nk = K.q(k); lastk = firstk + nk - 1;
     z(firstk:lastk) = qjmul(x(firstk:lastk), y(firstk:lastk));
     firstk = lastk + 1;
   end
 end
% ----------------------------------------
% SDP:  Z = (XY + YX)/2  if x full, or
% ----------------------------------------
 if length(x) == length(y)
   for k = 1:K.rsdpN                               % REAL symmetric PSD
     nk = K.s(k); lastk = firstk + nk*nk -1;
     Xk = mat(x(firstk:lastk),nk);
     Yk = mat(y(firstk:lastk),nk);
     z(firstk:lastk) = (Xk'*Yk+ Yk'*Xk)/2;
     firstk = lastk + 1;
   end
   for k = (1+K.rsdpN):length(K.s)                  % COMPLEX Hermitian PSD
     nk = K.s(k); imfirstk = firstk + nk*nk; lastk = imfirstk + nk*nk -1;
     Xk = mat(x(firstk:imfirstk-1) + sqrt(-1)*x(imfirstk:lastk),nk);
     Yk = mat(y(firstk:imfirstk-1) + sqrt(-1)*y(imfirstk:lastk),nk);
     Zk = (Xk'*Yk+ Yk'*Xk)/2;
     z(firstk:imfirstk-1) = real(Zk);
     z(imfirstk:lastk) = imag(Zk);
     firstk = lastk + 1;
   end
 else
% ----------------------------------------
% SDP:  z(i.j) = y(i,j)*(xi + xj)/2 if x spectral.
% ----------------------------------------
   for k = 1:K.rsdpN                               % REAL symmetric PSD
     nk = K.s(k);
     xk = x(nextx:nextx+nk-1);
     for j = 1:nk
       lastk = firstk + nk - 1;
       xkj = xk(j);
       z(firstk:lastk) = y(firstk:lastk).*(xkj+xk)/2;
       firstk = lastk + 1;
     end
     nextx = nextx + nk;
   end
   for k = (1+K.rsdpN):length(K.s)                  % COMPLEX Hermitian PSD
     nk = K.s(k);
     xk = x(nextx:nextx+nk-1);
     for j = 1:nk                                   % REAL part of y,z
       lastk = firstk + nk - 1;
       xkj = xk(j);
       z(firstk:lastk) = y(firstk:lastk).*(xkj+xk)/2;
       firstk = lastk + 1;
     end
     for j = 1:nk                                    % IMAGINARY part of y,z
       lastk = firstk + nk - 1;
       xkj = xk(j);
       z(firstk:lastk) = y(firstk:lastk).*(xkj+xk)/2;
       firstk = lastk + 1;
     end
     nextx = nextx + nk;
   end
 end
