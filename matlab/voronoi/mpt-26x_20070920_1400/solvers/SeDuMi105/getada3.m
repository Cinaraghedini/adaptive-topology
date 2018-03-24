%                            [ADA,absd] = getada3(ADA, A,Ajc1,Aord, udsqr,K)
% GETADA3  Compute ADA(i,j) = (D(d^2)*A.t(:,i))' *A.t(:,j),
%   and exploit sparsity as much as possible.
%   absd - length m output vector, containing
%     absd(i) = abs((D(d^2)*A.t(:,i))' *abs(A.t(:,i)).
%     Hence, diag(ADA)./absd gives a measure of cancelation (in [0,1]).
%
% SEE ALSO sedumi, getada1, getada2
% ******************** INTERNAL FUNCTION OF SEDUMI ********************
function [ADA,absd] = getada3(ADA, A,Ajc1,Aord, udsqr,K)
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
