function [Xid,yid]=hit_pwarx_format_data(u,y,na,nb);
%HIT_PWARX_FORMAT_DATA Build regressors and outputs for a MISO PWARX model of orders na and nb.
%
% -------------------------------------------------------------------------
% DESCRIPTION
% -------------------------------------------------------------------------
%
% [Xid,yid]=hit_pwarx_format_data(u,y,na,nb);
%
% Each mode of the PWARX system is assumed to be in the form
% y(k)=theta_i * [y(k-1) ... y(k-na) u'(k-1) ... u'(k-nb)]' + constant
%
% -------------------------------------------------------------------------
% INPUT
% -------------------------------------------------------------------------
% u(i,:) is the i-th input sample.
%
% y(i) is the i-th output sample.
%
% na, nb: orders of the PWARX model.
%
% -------------------------------------------------------------------------
% OUTPUT                                                                                                   
% -------------------------------------------------------------------------
% Xid: matrix containing the regressors. Each row is a datapoint.
%
% yid: column vector containing the output datapoints.

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

Xid=[];
yid=[];

% put outputs in column format
y=y(:);
N=length(y);
dimu=size(u,2);
Nu=size(u,1);
if Nu<=dimu
    fprintf('\n ======================================================== ');
    fprintf('\n DATA CHECK ERROR: insufficient number of input datapoints');
    fprintf('\n Check if the inputs are stored as ROWS of u. ');
    fprintf('\n ======================================================== \n');
    error('hit_pwarx_format_data:end','hit_pwarx_format_data ends here.')
end
if Nu>=N;
    fprintf('\n ======================================================== ');
    fprintf('\n DATA CHECK ERROR: the number of input samples is >=');
    fprintf('\n the number of output samples. ');
    fprintf('\n ======================================================== \n');
    error('hit_pwarx_format_data:end','hit_pwarx_format_data ends here.')
end
% mn: maximum n. of input/output samples in the past
mn=max(na,nb);
if mn<=0 | min(na,nb)<0
    fprintf('\n ======================================================== ');
    fprintf('\n DATA CHECK ERROR: the system order are both <=0');
    fprintf('\n or one system order is negative');
    fprintf('\n ======================================================== \n');
    error('hit_pwarx_format_data:end','hit_pwarx_format_data ends here.')
end
for k=mn+1:N
    reg=[];
    for j=1:na
        reg=[reg y(k-j)];
    end
    for j=1:nb
        reg=[reg u(k-j,:)];
    end
    Xid=[Xid;reg];
    yid=[yid ; y(k)];
end
