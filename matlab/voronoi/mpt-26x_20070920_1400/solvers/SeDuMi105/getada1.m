%                                 ADA = getada1(ADA, A,Ajc2,perm, d, blkstart)
% GETADA1  Compute ADA(i,j) = (D(d^2; LP,Lorentz)*A.t(:,i))' *A.t(:,j),
%   and exploit sparsity as much as possible.
%   Ajc2 points just beyond LP/Lorentz nonzeros for each column
%   blkstart = K.qblkstart partitions into Lorentz blocks.
%
%   IMPORTANT 1: only LP and sparse Lorentz part. PSD part ignored altogether.
%   For Lorentz, it uses only det(dk) * ai[k]'*aj[k].
%   IMPORTANT 2: Computes ADA only on triu(ADA(Aord.lqperm,Aord.lqperm)).
%     Remaining entries are set to 0. (CAUTION: sparse(ADA) will therefore
%     destroy the sparsity structure !).
%
% SEE ALSO sedumi, getada2, getada3
% ******************** INTERNAL FUNCTION OF SEDUMI ********************
function ADA = getada1(ADA, A,Ajc2,perm, d, blkstart)
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

error('At OS prompt, type "make" to create SeDuMi mex-files.')
