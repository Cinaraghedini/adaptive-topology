%                                                y = Amul(At,dense,x,transp)
% AMUL  Computes A*x (transp=0) or A'*x (transp=1), taking care of dense.A.
%
% SEE ALSO sedumi
% **********  INTERNAL FUNCTION OF SEDUMI **********
function y = Amul(At,dense,x,transp)
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

if nargin < 4
  transp = 0;
end
if transp == 0
  y = At'*x;                       %y(m) and x(N).
else
  y = full(At*x);        % y(N) and x(m)
end
if length(dense.cols) > 0
  if transp == 0
    y = y + dense.A*x(dense.cols);
  else
    y = full(At*x);                
    y(dense.cols) = dense.A'*x;
  end
end
