% cpx = complexblk(z,K)
% **********  INTERNAL FUNCTION OF SEDUMI **********
function cpx = complexblk(z,K)

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

% ----------------------------------------
% LP: should be real
% ----------------------------------------
if ~isreal(z(1:K.l))
  error('Complex data for nonnegative real variables encountered')
end
% ----------------------------------------
% LORENTZ: find complex blocks
% ----------------------------------------
cpx.q = zeros(length(K.q),1);
firstk = K.l + 1;
for k = 1:length(K.q)
  nk = K.q(k);
  if ~isreal(z(firstk))
    error('Complex data in Lorentz real scalar part encountered')
  end
  cpx.q(k) = ~isreal(z(firstk+1:firstk+nk-1));
  firstk = firstk + nk;
end
cpx.q = find(cpx.q);
% ----------------------------------------
% PSD: 
% ----------------------------------------
cpx.s = zeros(length(K.s),1);
for k = 1:length(K.s)
  nksqr = K.s(k)^2;
  cpx.s(k) = ~isreal(z(firstk+1:firstk+nksqr-1));
  firstk = firstk + nksqr;
end
cpx.s = find(cpx.s);
% ----------------------------------------
% Let cpx.dim = number of imaginary-valued variables
% ----------------------------------------
cpx.dim = (sum(K.q(cpx.q)) - length(cpx.q)) + sum(K.s(cpx.s).^2);
