%                                      [At,b,c,K,prep] = pretransfo(At,b,c,K)
% PRETRANSFO  Checks data and then transforms into internal SeDuMi format.
%
% SEE ALSO sedumi
% **********  INTERNAL FUNCTION OF SEDUMI **********
function [At,b,c,K,prep] = pretransfo(At,b,c,K)

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

% ----------------------------------------
% Make sure that all fields exist in K-structure
% ----------------------------------------
  if ~isfield(K,'f')                       % K.f
    K.f = 0;
  elseif isempty(K.f)
    K.f = 0;
  elseif (K.f~= floor(K.f)) | (K.f < 0)
    error('K.f should be nonnegative integer')
  end
  if ~isfield(K,'l')                        % K.l
    K.l = 0;
  elseif isempty(K.l)
    K.l = 0;
  elseif (K.l~= floor(K.l)) | (K.l < 0)
    error('K.l should be nonnegative integer')
  end
  if ~isfield(K,'q')                       % K.q
    K.q = [];
  elseif sum(K.q) == 0
    K.q = [];
  elseif ~isempty(K.q)
    if (K.q~= floor(K.q)) | min(K.q) < 2
      error('K.q should contain only integers bigger than 1')
    end
    if size(K.q,1) > 1
      K.q = K.q';
    end
  end
  if ~isfield(K,'r')                       % K.r
    K.r = [];
  elseif sum(K.r) == 0
    K.r = [];
  elseif ~isempty(K.r)
    if (K.r~= floor(K.r)) | min(K.r) < 3
      error('K.r should contain only integers bigger than 2')
    end
    if size(K.r,1) > 1
      K.r = K.r';
    end
  end
  if ~isfield(K,'s')                       % K.s
    K.s = [];
  elseif sum(K.s) == 0
    K.s = [];
  elseif ~isempty(K.s)
    if (K.s~= floor(K.s)) | min(K.s < 1)
      error('K.s should contain only positive integers')
    end
    if size(K.s,1) > 1
      K.s = K.s';
    end
  end
% ------------------------------------------------------------
% Validity check: Info on complex interpretation of variables
% ------------------------------------------------------------
  if ~isfield(K,'ycomplex')                        % K.ycomplex
    K.ycomplex = [];
  elseif ~isempty(K.ycomplex)
    if (K.ycomplex~= floor(K.ycomplex)) | min(K.ycomplex) < 1 | ...
        min(size(K.ycomplex)) > 1
      error('K.ycomplex should be a list containing only positive integers')
    end
    if size(K.ycomplex,1) > 1
      K.ycomplex = K.ycomplex';
    end
    if max(K.ycomplex) > length(b)
      error('K.ycomplex out of range')
    end
  end
  if ~isfield(K,'xcomplex')                        % K.xcomplex
    K.xcomplex = [];
  elseif ~isempty(K.xcomplex)
    if (K.xcomplex~= floor(K.xcomplex)) | min(K.xcomplex) < 1 | ...
        min(size(K.xcomplex)) > 1
      error('K.xcomplex should be a list containing only positive integers')
    end
    if size(K.xcomplex,1) > 1
      K.xcomplex = K.xcomplex';
    end
    if max(K.xcomplex) > K.f+K.l+sum(K.q)+sum(K.r)
      error('K.xcomplex out of range')
    end
  end
  if ~isfield(K,'scomplex')                        % K.scomplex
    K.scomplex = [];
  elseif ~isempty(K.scomplex)
    if (K.scomplex~= floor(K.scomplex)) | min(K.scomplex) < 1 | ...
        min(size(K.scomplex)) > 1
      error('K.scomplex should be a list containing only positive integers')
    end
    if size(K.scomplex,1) > 1
      K.scomplex = K.scomplex';
    end
    if max(K.scomplex) > length(K.s)
      error('K.scomplex out of range')
    end
  end
% ----------------------------------------
% Check size of At,b,c (w.r.t. K)
% Let m = #eq-constraints, N = #variables (before transformations)
% ----------------------------------------
  if (min(size(b)) > 1) | (min(size(c)) > 1)
    error('Parameters b and c must be vectors')
  end
  m = min(size(At));
  N = K.f + K.l + sum(K.q) + sum(K.r) + sum((K.s).^2);
  if nnz(b) == 0
    b = sparse(m,1);
  elseif length(b) ~= m
    error('Size b mismatch')
  end
  if m ~= size(b,1)                    % make b an m x 1 vector
    b = b';
  end
  if nnz(c) == 0
    c = sparse(N,1);
  elseif length(c) ~= N
    error('Size c mismatch')
  end
  if N ~= size(c,1)                     % make c an N x 1 vector
    c = c';
  end
  if N <= m
    error('Should have length(c) > length(b) in any sensible model')
  end
  if size(At,2) ~= m 
    if m == size(At,1)
      At = At';        %user gave A instead of At.
    else
      error('Size At mismatches b.')
    end
  end
  if size(At,1) ~= N;
    error('Size At mismatches cone K.')
  end
% ------------------------------------------------------------
% Check for NaN and Inf
% ------------------------------------------------------------
 if any(any(isnan(At))) | any(isnan(b)) | any(isnan(c))
   error('A,b,c data contains NaN values')
 end
 if any(any(isinf(At))) | any(isinf(b)) | any(isinf(c))
   error('A,b,c data contains Inf values')
 end
% ----------------------------------------
% First add constraints IM(Ax)=IM(b) where i in K.ycomplex, i.e.
% Make "imag(Ax) = imag(b)" and "imag(y)" explicit.
% ----------------------------------------
  if ~isempty(K.ycomplex)
    b = [real(b);...
         imag(b(K.ycomplex))];      % RE b'*y = [RE b; IM b]'*[RE y; IM y].
    At = [At, sqrt(-1) * At(:,K.ycomplex)];
  else
    b = real(b);
  end
  m = length(b);
% ----------------------------------------
% Transform any complex data into internal format,
% which uses only MATLAB's real representation.
% ----------------------------------------
  cpx = whichcpx(K);
  cpx.s = K.scomplex;
  cpx.dim = length(K.xcomplex) + sum(K.s(K.scomplex).^2);
  At = makereal(sparse(At),K,cpx);
  c = makereal(sparse(c),K,cpx);
  K.f = K.f + length(cpx.f);          % Update cone K structure wrt complex
  K.q = vec(K.q) + vec(cpx.q);
  K.r = vec(K.r) + vec(cpx.r);
  sperm = ones(length(K.s),1);
  sperm(cpx.s) = 0;
  sperm = find(sperm);
  K.rsdpN = length(sperm);          % #real sym PSD blocks
  K.s = K.s([vec(sperm); vec(cpx.s)]);
  prep.cpx = cpx;
% ------------------------------------------------------------
% Transform F-cone into L-cone
% ------------------------------------------------------------
  prep.Kf = K.f;
  At = [-At(1:K.f,:);At];
  c = [-c(1:K.f);c];
  K.l = K.l + 2 * K.f;
  K.f = 0;
% ----------------------------------------
% Transform R-cones (rotated Lorentz) into Q-cones (standard Lorentz).
% ----------------------------------------
  prep.lenKq = length(K.q);
  if ~isempty(K.r)
    c = rotlorentz(c,K);
    At = rotlorentz(At,K);
    K.q = [K.q; K.r];
    K.r = [];
  end
% ----------------------------------------
% Now that complex Lorentz cones are transformed to real,
% we need K.q(k) >= 3 for all k.
% ----------------------------------------
  if min(K.q) < 3
    error('The real dimension of each Lorentz cone should be at least 3.')
  end
% --------------------------------------------------
% Split Lorentz blocks in trace-block + norm-bound-blocks
% --------------------------------------------------
  At = qreshape(At,0,K);
  c = qreshape(c,0,K);
% ----------------------------------------
% Correct A s.t. At has tril-blocks (for PSD),
% ----------------------------------------
  At = vectril(At,K);
  c = vectril(sparse(c),K);
% ----------------------------------------
% Create artificial (x0,z0) variable for
% self-dual model
% ----------------------------------------
  c = [0;c];            %does not affect sparse/dense status of c.
  At = [sparse(1,m); At];
  K.l = K.l + 1;    % add (x0,z)
% ----------------------------------------
% Now K has field K.{l,q,s}
% Make more detailed description of cone K:
% Let K.blkstart(i):K.blkstart(i+1)-1 give the index range of block i.
% Compute maxN=max(order), len=sum(order) for LORENTZ, real PSD and herm PSD.
% yields: K.{l,q,s, rsdpN,blkstart,    rLen,hLen, qMaxn,rMaxn,hMaxn}
% ----------------------------------------
  K.rsdpN = length(K.s) - length(prep.cpx.s);    % # real symmetric PSD vars
  K = statsK(K);
  K.mainblks = K.blkstart(cumsum([1 1 length(K.q)]));
  K.qblkstart = K.blkstart(2:2+length(K.q));  % Also include blkend
  K.sblkstart = K.blkstart(2+length(K.q):end);
  K.lq = K.mainblks(end)-1;
  K.N = length(c);
