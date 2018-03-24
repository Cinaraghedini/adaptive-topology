% z = invjmulK(x,y, K,f)
% solves x jmul z = y, with x in spectral space and Lorentz frame f.
% **********  INTERNAL FUNCTION OF SEDUMI **********

function z = invjmulK(x,y, K,f)

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

%error('At OS prompt, type "make" to create SeDuMi mex-files.')

 z = zeros(length(y),1);
% ----------------------------------------
% LP: z = y./x
% ----------------------------------------
 z(1:K.l) = y(1:K.l)./x(1:K.l);
 firstk = K.l + 1;
 nextx = K.l+1; nextf = 1;
% ----------------------------------------
% LORENTZ:
% tr z = inv(x)'*y
% z_2 = (2*y_2 - (tr z)*x_2) / tr x.
% ----------------------------------------
 for k=1:length(K.q)
   nk = K.q(k); lastk = firstk + nk - 1;
   fy2 = f(nextf:nextf+nk-2)' * y(firstk+1:lastk);
   x1r2 = (x(nextx) + x(nextx+1)) / 2;
   difxr2 = (x(nextx+1) - x(nextx)) / sqrt(2);
   z1 = (x1r2 * y(firstk) - difxr2 * fy2) / (x(nextx)*x(nextx+1));
   z2 = y(firstk+1:lastk)/x1r2 - (z1 * difxr2 / x1r2) * f(nextf:nextf+nk-2);
   z(firstk) = z1;
   z(firstk+1:lastk) = z2;
   firstk = lastk + 1; nextf = nextf + nk-1; nextx = nextx + 2;
 end
% ----------------------------------------
% real PSD:  z(i.j) = 2*y(i,j)/(xi + xj)   (Z and Y are matrices, x is vector)
% ----------------------------------------
for k = 1:K.rsdpN
  nk = K.s(k);
  xk = x(nextx:nextx+nk-1);
  for j = 1:nk
    lastk = firstk + nk - 1;
    xkj = xk(j);
    z(firstk:lastk) = 2*y(firstk:lastk)./(xkj+xk);
    firstk = lastk + 1;
  end
  nextx = nextx + nk;
end
% ----------------------------------------
% Hermitian PSD:  z(i.j) = 2*y(i,j)/(xi + xj)
%  Z and Y are Hermitian matrices, x is real vector
% ----------------------------------------
for k = (1+K.rsdpN):length(K.s)
  nk = K.s(k);
  xk = x(nextx:nextx+nk-1);
  for j = 1:nk                         % 1: real part
    lastk = firstk + nk - 1;
    xkj = xk(j);
    z(firstk:lastk) = 2*y(firstk:lastk)./(xkj+xk);
    firstk = lastk + 1;
  end
  for j = 1:nk                         % 2: imaginary part
    lastk = firstk + nk - 1;
    xkj = xk(j);
    z(firstk:lastk) = 2*y(firstk:lastk)./(xkj+xk);
    firstk = lastk + 1;
  end
  nextx = nextx + nk;
end
