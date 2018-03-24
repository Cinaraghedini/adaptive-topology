% GETPROBLEM  Script for locating, reading and convert to SeDuMi, a problem
%   with the name specified by the variable `pname'. You have to define the
%   variable 'pname' before using this script.
%
% Usage: pname = 'problemname'; getproblem
%   Problem returned in (At,b,c,K).
%
% IMPORTANT
%   1. If you have an environment variable named SDPPACK, containing
%   the directory name of SDPPACK, then this searches in SDPPACK/problems.
%   2. If you have an environment variable named SDPLIB, then it searches
%   for problems in the directory SDPLIB/. These problems should be in
%   sparse SDPA format, with extension ".dat-s".
%   3. If LIPSOL is installed as a MATLAB Toolbox, then LIPSOL's  `findprob'
%   is used for searching MPS-type models.
%   4. Files with extension ".gz" are also recognized, so that you can
%   store your problem libraries more efficiently.
%
% SEE ALSO fromsdpa, frompack, prelp, feasreal, feascpx, sedumi,
%   postprocess (LIPSOL).

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

global PROBDIR

% ----------------------------------------
% status = {0 don't know yet, 1 SDPPack, 2 MPS, 3 SDPA}
% ----------------------------------------
status = 0;
PROBDIR = str2mat([]);
% ----------------------------------------
% Try to get problem in SDPA format, from SDPLIB
% ----------------------------------------
if status == 0
  SDPAPATH = getenv('SDPLIB');
  if ~isempty(SDPAPATH)
    MATNAME = [SDPAPATH '/' pname '.dat-s'];
    if exist(MATNAME) 
      [At,b,c,K] = fromsdpa(MATNAME); status = 3;
    elseif exist([MATNAME '.gz'])
      unix(['gunzip -c ' MATNAME '.gz > /tmp/default.mat']);
      [At,b,c,K] = fromsdpa('/tmp/default.mat'); status = 3;
    end
  end
  clear SDPAPATH
end
% ----------------------------------------
% Try to get problem in SDPPACK format
% ----------------------------------------
if status == 0
  PACKPATH = getenv('SDPPACK');
  if ~isempty(PACKPATH)
    MATNAME = [PACKPATH '/problems/' pname '.mat'];
    if exist(MATNAME) 
      unix(['/bin/cp ' MATNAME ' /tmp/default.mat']); status = 1;
    elseif exist([MATNAME '.gz'])
      unix(['gunzip -c ' MATNAME '.gz > /tmp/default.mat']); status = 1;
    end
  end
  clear PACKPATH
end
% ----------------------------------------
% Try to get problem in LIPSOL format
% ----------------------------------------
if(status == 0)
  if(exist('findprob') > 1)
    if (findprob([],pname) == 1)
      status = 2;
    end
  end
end
% ----------------------------------------
% Move the problem into [At,b,c,K] structure
% ----------------------------------------
if(status==0)
  error('Problem not found.  Check spelling?')
elseif(status==1)
  load /tmp/default
% ----------------------------------------
% In SDPPACK, 0s are sometimes used  as emptyset and vice versa.
% ----------------------------------------
  if isempty(blk.l)
    blk.l = 0;
  end
  if ~isempty(blk.q)
    if(blk.q == 0)
      blk.q = [];
    end
  end
  if ~isempty(blk.s)
    if(blk.s == 0)
      blk.s = [];
    end
  end
  [At,c] = frompack(A,b,C,blk);
  if isempty(blk.l)
    blk.l = 0;
  end
  K = blk;
% ----------------------------------------
% LP problem: preprocess using LIPSOL.
% ----------------------------------------
elseif(status==2)
  [At,b,c,nx,lbs] = prelp([]);
  At = At';
  clear K;
  K.l = length(c);
end
clear status
% clear global PROBDIR   %would disable LIPSOL.