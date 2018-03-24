% [L,perm,xsuper,split,tmpsiz] = symfctmex(X, perm, cachsz)
%   Computes sparse symbolic factor L, updated permutation PERM,
%   super-node partition XSUPER, and a splitting of supernodes
%   (SPLIT) to optimize use of the computer cache (assuming
%   CACHSZ*1024 byte available).   TMPSIZ is the amount of floating
%   point working storage that has to be allocated within blkfctmex.
%
%   Invokes ORNL block Cholesky library (Fortran).
% **********  INTERNAL FUNCTION OF CHOLTOOL  **********

function [L,perm,xsuper,split,tmpsiz] = symfctmex(adjncy, perm, cachsz)
% THE M-FILE VERSION OF THIS FUNCTION IS HERE ONLY AS ILLUSTRATION.
% SEE THE F-SOURCE FOR THE MEX-VERSION.

 %  
 %   This file is part of CholTool 1.00
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

error('At OS prompt, type "make" to create cholTool mex-files.')



