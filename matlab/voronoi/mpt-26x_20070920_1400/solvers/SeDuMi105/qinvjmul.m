%                                          y = qinvjmul(labx,frmx,b,K)
% QINVJMUL  Inverse of Jordan multiply for Lorentz blocks
%
% SEE ALSO sedumi
% **********  INTERNAL FUNCTION OF SEDUMI **********
function y = qinvjmul(labx,frmx,b,K)
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
 lorN = length(K.q);
 if lorN == 0
   y = zeros(0,1);
   return
 end
 if length(labx) > 2*lorN
   labx = labx(K.l+1:K.l+2*lorN);
 end
 detx = labx(1:lorN) .* labx(lorN+1:end);
 x = qframeit(labx,frmx,K);
 ix = K.mainblks;
 if length(b) == ix(3)-ix(1);   % lorentz only ?
   ix = (1-ix(1)) + ix;
 end
% ------------------------------------------------------------
% Let y1(k) = xk'Jbk/(sqrt2*detxk)
% ------------------------------------------------------------
 y1 = x(1:lorN).*b(ix(1):ix(2)-1) - ddot(x(lorN+1:end),b,K.qblkstart);
 y1 = y1./(sqrt(2)*detx);
% ------------------------------------------------------------
% Let y2[k] = (sqrt2/x1)*b2[k] - (y1/x1) * x2[k]
% ------------------------------------------------------------
 y = [y1; qblkmul(sqrt(2)./x(1:lorN),b,K.qblkstart)...
       - qblkmul(y1./x(1:lorN),x(lorN+1:end),K.qblkstart)];
