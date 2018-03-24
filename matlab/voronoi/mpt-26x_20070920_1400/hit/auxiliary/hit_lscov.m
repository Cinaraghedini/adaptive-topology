function [x,stdx,mse,S] =hit_lscov(A,b,V)
%HIT_LSCOV Call to MatLab lscov.m in two different ways if Matlab <14 or >=14
%
% -------------------------------------------------------------------------
% DESCRIPTION
% -------------------------------------------------------------------------
% [x,stdx,mse,S] =hit_lscov(A,b,V)
%
% -------------------------------------------------------------------------
% INPUT
% -------------------------------------------------------------------------
% The same as lscov: type 'help lscov'. The only difference is that V is
% a vetor of weights (corresponding to the covariance matrix diag(V)).
%
% -------------------------------------------------------------------------
% OUTPUT                                                                                                   
% -------------------------------------------------------------------------
% The same as lscov in Matlab >=14, i.e.
%
% estimate: X = inv(A'*inv(V)*A)*A'*inv(V)*b
% Mean Squared Error: MSE = b'*(inv(V) - inv(V)*A*inv(A'*inv(V)*A)*A'*inv(V))*b./(M-N)
% estimate covariance: S = inv(A'*inv(V)*A)*MSE
% estimate standard deviation: STDX = sqrt(diag(S))

% Copyright is with the following author:
%
% (C) 2005 Giancarlo Ferrari Trecate,
%         giancarlo.ferrari@unipv.it
% -------------------------------------------------------------------------
% Legal note:
%     This program is free software; you can redistribute it and/or
%     modify it under the terms of the GNU General Public
%     License as published by the Free Software Foundation; either
%     version 2.1 of the License, or (at your option) any later version.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%     General Public License for more details.
%
%     You should have received a copy of the GNU General Public
%     License along with this library; if not, write to the
%     Free Software Foundation, Inc.,
%     59 Temple Place, Suite 330,
%     Boston, MA  02111-1307  USA
%
%
% -------------------------------------------------------------------------

if str2num(version('-release'))>=17
    [x,stdx,mse,S] =lscov(A,b,V);
else
    [x,stdx] = lscov(A,b,diag(V));
    % compute the inverse of the DIAGONAL matrix V
    IV=diag(1./V);
    IAIVA=inv(A'*IV*A);
    M=size(A,1);
    N=size(A,2);
   mse = b'*(IV - IV*A*IAIVA*A'*IV)*b./(M-N);
   S = IAIVA*mse;
end