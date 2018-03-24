function [best_s,best_c, mse_m]=hit_estimate_cs(Xid,yid,Xv,yv,c_test,s_test)
%HIT_ESTIMATE_CS Estimate the best mode number and/or the best size of LDs within a given range 
%
% -------------------------------------------------------------------------
% DESCRIPTION
% -------------------------------------------------------------------------
%
% [idmodes,F,xi,LDs,inliers]=hit_regression(Xid,yid)
% 
% -------------------------------------------------------------------------
% INPUT
% -------------------------------------------------------------------------
% Xid: matrix containing the regressors. Each row is a datapoint.
%
% yid: column vector containing the output datapoints.
%
% Xv: matrix containing validation regressors. Each row is a datapoint.
%
% yv: column vector containing the validation output datapoints.
%
% c_test: vector of integers collecting the values of idpar.c (the LD
% size) to try.
%
% s_test: vector of integers collecting the values of idpar.s (the number
% of modes) to try.
%
% -------------------------------------------------------------------------
% OUTPUT                                                                                                   
% -------------------------------------------------------------------------
% best_c, best_s: the joint combination of LD size and mode number that
% minimizes the MSE on validation data.
%
% mse_m(i,j): MSE on validation data when using LD size equal to c_test(i)
% and number of modes equal to s_test(i).

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

global idpar plotpar ;
if ~isstruct(idpar) | ~isstruct(plotpar),
    hit_error;
end

% save some fileds of idpar and plotpar
% to be restored at the end
c_saved=idpar.c;
s_saved=idpar.s;
plot_class_xi_points_yn_saved=plotpar.plot_class_xi_points_yn;
plot_xi_points_yn_saved=plotpar.plot_xi_points_yn;
patt_rec_algo_saved=idpar.patt_rec_algo;

% Options
plotpar.plot_xi_points_yn='N';
plotpar.plot_class_xi_points_yn='N';

% Use a fast pattern recognition algorithm
idpar.patt_rec_algo='PSVC';

% Main loop
nc=length(c_test);
ns=length(s_test);
mse_m=zeros(nc,ns);
ndatav=length(yv);
for i=1:nc
    idpar.c=c_test(i);
    for j=1:ns
        idpar.s=s_test(j);
        fprintf('\n ======================================================== ');
        fprintf('\n Building model with c=%d and s=%d modes',c_test(i),s_test(j));
        fprintf('\n ======================================================== \n');
        [idmodes{i,j},F_trash,xi_trash,LDs_trash,inliers_trash]=hit_regression(Xid,yid);
        % compute the mean square error 
        mse_m(i,j)=hit_mse(Xv,yv,idmodes{i,j});
    end
end

% find best model according to the mse on validation data

best_mse=min(min(mse_m));
%find the position of best_ssr in ssr_m
[p,q]=find(mse_m==ones(nc,ns)*best_mse);
%restore the original fields of idpar and plotpar
idpar.c=c_saved;
idpar.s=s_saved;
plotpar.plot_class_xi_points_yn=plot_class_xi_points_yn_saved;
plotpar.plot_xi_points_yn=plot_xi_points_yn_saved;
idpar.patt_rec_algo=patt_rec_algo_saved;

% outputs

best_c=c_test(p);
best_s=idmodes{p,q}.s;
best_model=idmodes{p,q};
fprintf('\n ======================================================== ');
fprintf('\n Best parameters: c=%d and s=%d ',best_c,best_s);
fprintf('\n ======================================================== \n');
