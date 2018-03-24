% Y = MAT(x,n)   or   Y = MAT(x)     (the 2nd argument is optional)
%   Given a vector of length n^2, this produces the n x n matrix
%   Y such that x = vec(Y).  In other words, x contains the columns of the
%   matrix Y, stacked below each other.
%
% SEE ALSO vec.

function X = mat(x,n)

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

 if nargin < 2
   n = floor(sqrt(length(x)));
   if (n*n) ~= length(x)
     error('Argument X has to be a square matrix')
   end
 end
 X = reshape(x,n,n);
%if issparse(x)
%  X = sparse(1+rem(find(x)-1,n), 1+fix((find(x)-1)/n), nonzeros(x),n,n);
%else
%  X = zeros(n,n);
%  X(:) = full(x);
%  if issparse(x)
%    X = sparse(X);
%  end
%end
