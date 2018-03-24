% c = rotlorentz(c,K)
% Rotates vectors from Qcone to Rcone or from Rcone into Qcone.
% **********  INTERNAL FUNCTION OF SEDUMI **********
function c = rotlorentz(c,K)

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

firstk = K.l + sum(K.q) + 1;
M = [1 1; 1 -1] / sqrt(2);
for k = 1:length(K.r)
  c1 = c(firstk); c2 = c(firstk+1);
  c(firstk:firstk+1,:) = M*c(firstk:firstk+1,:);
  firstk = firstk + K.r(k);
end
