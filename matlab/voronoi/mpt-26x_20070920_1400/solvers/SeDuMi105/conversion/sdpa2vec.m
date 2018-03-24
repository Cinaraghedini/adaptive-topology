% x = sdpavec(E,K)
% Takes an SDPA type sparse data description E, i.e.
%   E(1,:) = block, E(2,:) = row, E(3,:) = column, E(4,:) = entry,
% and transforms it into a "long" vector, with vectorized matrices for
% each block stacked under each other. The size of each matrix block
% is given in the field K.s.
% **********  INTERNAL FUNCTION OF SEDUMI **********

function x = sdpa2vec(E,K,invperm)

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
% Split E into blocks E[1], E[2], ..., E[length(K.s)]
% ------------------------------------------------------------
 [xir,xjc] = sdpasplit(E(1,:),length(K.s));
 blkstart = cumsum([K.l+1; K.s.^2]);
 x = sparse([],[],[],K.l+sum(K.s.^2),1,2*size(E,2));
 for knz = 1:length(xir)
   permk = xir(knz);
   k = invperm(permk);
   if k > K.l
     k = k - K.l;                       % matrix
     nk = K.s(k);
     Ek = E(2:4,xjc(knz):xjc(knz+1)-1);
     Xk = sparse(Ek(1,:),Ek(2,:),Ek(3,:),nk,nk);
     Xk = Xk + triu(Xk,1)';
     x(blkstart(k):blkstart(k+1)-1) = vec(Xk);
   else
     x(k) = E(4,xjc(knz));      % Scalar
   end
 end
