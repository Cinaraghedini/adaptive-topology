% PRELP  Loads and preprocesses LP from an MPS file.
%
% > [A,b,c,lenx,lbounds] = PRELP('problemname')
%    The above command results in an LP in standard form,
%    - Instead of specifying the problemname, you can also use PRELP([]), to
%    get the problem from the file /tmp/default.mat.
%    - Also, you may type PRELP without any input arguments, and get prompted
%    for a name.
%
% MINIMIZE  c'*x SUCH THAT  A*x = b AND  x>= 0.
%
%   So, you can solve it with SeDuMi:
% > [x,y,info] = SEDUMI(A,b,c);
%
% After solving, post-process it with
% > [x,objp] = POSTPROCESS(x(1:lenx),lbounds).
%
% REMARK  x(lenx+1:length(x)) will contain upper-bound slacks.
%
% IMPORTANT works only if LIPSOL is installed on your system.
%
% SEE ALSO sedumi, getproblem, postprocess (LIPSOL), frompack AND lipsol.

function [A,b,c,lenx,lbounds] = prelp(pname)
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


global OUTFID
global Ubounds_exist

if ~exist('loadata') | ~exist('preprocess')
  error('To use PRELP, you need to have LIPSOL installed.')
end

%--------------------------------------------------
% LOAD LP PROBLEM INTO MEMORY
%--------------------------------------------------
if (nargin == 0) pname = input('Enter problem name: ','s'); end

t0 = cputime;
[A,b,c,lbounds,ubounds,BIG,NAME] = loadata(pname);
times(1) = cputime - t0;

%--------------------------------------------------
% PREPROCESS LP PROBLEM
% NB: Y.Zhang's preprocess returns lbounds for post-
% processing; the pre-processed problem has x>=0.
%--------------------------------------------------
t0 = cputime;
[A,b,c,lbounds,ubounds,FEASIBLE] = ...
   preprocess(A,b,c,lbounds,ubounds,BIG);
if ~FEASIBLE
   fprintf('\n'); if isempty(OUTFID) return; end;
   msginf =  'Infeasibility detected in preprocessing';
   fprintf(OUTFID, [pname '  0   ' msginf '\n']);
   return;
end;
%[A,b,c,ubounds] = scaling(A,b,c,ubounds);
%--------------------------------------------------
% INSERT UBOUND- CONSTRAINTS IN THE A-MATRIX
%--------------------------------------------------
b = full(b); c = full(c);
[m,lenx] = size(A);
if Ubounds_exist
  nub = nnz(ubounds);
  A= [ A sparse(m,nub); sparse(1:nub,find(ubounds),1,nub,lenx) speye(nub) ];
  b = [b; nonzeros(ubounds)];
  c = [c; zeros(nub,1)];
else
  ubounds = [];
end
%--------------------------------------------------
% LOCATE DENSE COLUMNS
%--------------------------------------------------
%checkdense(A);
times(2) = cputime - t0;
