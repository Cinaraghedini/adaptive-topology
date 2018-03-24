function [x,y,info] = runsp(At,b,c,sdpL);
% [X,Y,INFO] = RUNSP(At,b,c,sdpL)
%    Runs Vandenberghe & Boyd's code "SP" (written in C)
%    for SDP min{ c'*x | At'*x=b, x in K} with K.s = sdpL, K.l=[], K.q=[].
%    See "sedumi" for an explanation of the K and INFO-structures.
% 
% IMPORTANT: this function assumes that the directory of SP
%    exists in your search path.
%
% SEE ALSO sedumi, frompack, fromsdpa, getproblem.

% Jos F. Sturm, Canada, 1998. Adapted from "lin_prog.m" by Vandenberghe
% and Boyd.

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


F = [full(c) full(At)];
b = full(b);
M = 1e3 * max(sum(abs(F)));
info.pinf = 0; info.dinf = 0;
info.failure = 0;

% JFS: VB advices 10 <= nu <= 50, so I take the middle:

nu = 30.0;
abstol = 1e-9;  reltol = 1e-7; maxiters = 100;

%
% phase 1
% 

disp(' '); disp(' PHASE 1.');

cputime0 = cputime;
[x0,Z0,z0,ul,iters1] = ...
   phase1(F,sdpL, M, nu, abstol, maxiters);

if iters1 == maxiters,
   info.failure = 1;
   x = []; y = [];
   return;
end;

if (ul(1) > 0.0) 
   info.dinf = 1;               % The dual (min b'*x | c + A*x >= 0) may be
   x=[]; y=[];                  % infeasible
   return;
end;


%
% phase 2
% 

disp(' '); disp(' PHASE 2.');
M = max(M, 1e3*sum(c+At*x0));    % M must be greater than e'(b-A*x0) 
                                % for bigM.m 
[x,z,w,ul,iters2] = ...
   bigM(F, sdpL, b, x0, M, nu, abstol, reltol, 0.0, maxiters);

info.cpusec = cputime-cputime0;
info.iter = iters1 + iters2;

if iters2 == maxiters
   info.failure = 1;
   x=[]; y=[];
   return;
end;

if sum(F*[1;x]) > 0.9*M                     % This is a rather strange rule:
   info.pinf = 1;                   % If the dual slack is big, then it could
   x=[]; y=[];                      % be unbounded, and so primal infeasible.
   return;
end;

% ------------------------------------------------------------
% SP's "x" is actually -y
% SP's "z" is actually x 
% ------------------------------------------------------------
y = -x;
x = z;
