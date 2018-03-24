% SPARBWSLV Solves block sparse upper-triangular system.
%    y = sparbwslv(L,b) yields the same result as
%              y(L.perm,:) = L.L'\b
%    However, SPARBWSLV is faster than the built-in operator "\",
%    because it uses dense linear algebra and loop-unrolling on
%    supernodes.
%
%    Typical use, with X sparse m x m positive definite and b is m x n:
%            L = sparchol(symbchol(X),X);
%            L.d(L.dep) = inf;
%            y = sparbwslv(L,sparfwslv(L,b) ./ L.d);
%    Then y solves X*y=b.
%
% SEE ALSO symbchol, sparchol, sparfwslv, \, /.
function y = sparbwslv(L,b)
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

%% y2 = bwslvp(L.U22,L.ph2perm,b(L.ph2nodes));
%% y = b - L.U12 * y2;
%% y(L.ph2nodes) = y2;
 y = bwblkslv(L,b);
