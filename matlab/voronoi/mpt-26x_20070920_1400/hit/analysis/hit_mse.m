function [mse,mse_mode]=hit_mse(X,y,idmodes)
%HIT_MSE Estimate the Mean Squared Error of a PWA function on a given dataset
% -------------------------------------------------------------------------
% DESCRIPTION
% -------------------------------------------------------------------------
%
% [mse,mse_mode]=hit_mse(X,y,idmodes)
% 
% -------------------------------------------------------------------------
% INPUT
% -------------------------------------------------------------------------
% X: matrix containing the regressors. 
%
% y: column vector containing the output datapoints.
%
% idmodes: structure containing all information on the PWA function.
% Type 'help hit_regression' for a description of the fields. 
%
% -------------------------------------------------------------------------
% OUTPUT                                                                                                   
% -------------------------------------------------------------------------
%
% mse: mean squared error on the given data.
%
% mse_mode(i): mean squared error of the i-th mode on the given data. It is
% an optional output.

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

ndata=size(X,1);
y=y(:);
yp=zeros(ndata,1);
yr=yp;
for i=1:ndata
    [yp(i),yr(i)]=hit_pwa(idmodes.par,idmodes.regions,X(i,:)); % add noise to the data
    if isnan(yp(i))
        fprintf('\n ======================================================== ');
        fprintf('\n WARNING: the regressor [')
        fprintf('%f ',X(i,:))
        fprintf('] used for computing ssr');
        fprintf('\n falls outside all mode regions.');
        fprintf('\n The SSR for this point is set to zero by default');
        fprintf('\n ======================================================== \n');
        yp(i)=y(i);
    end
end

mse=(yp-y)'*(yp-y)/ndata;
if nargout>1
    for i=1:idmodes.s
        % indexes of points belonging to the i-th mode
        indpoints=find(yr==i);
        npoints=length(indpoints);
        mse_mode(i)=(yp(indpoints)-y(indpoints))'*(yp(indpoints)-y(indpoints))/npoints;
    end
end
    