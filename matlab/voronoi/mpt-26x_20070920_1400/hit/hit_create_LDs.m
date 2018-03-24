function LDs=hit_create_LDs(Xid,yid,Weight_primal)
%HIT_CREATE_LDS Create Local Datasets (LDs)
%
% -------------------------------------------------------------------------
% DESCRIPTION
% -------------------------------------------------------------------------
% LDs=hit_create_LDs(Xid,yid,Weight_primal)
%
% -------------------------------------------------------------------------
% INPUT
% -------------------------------------------------------------------------
%  Xid: matrix containing the datapoints in the regressor set.
%       Each row is a datapoint.
%
%  yid: column vector containing the output datapoints.
%
% -------------------------------------------------------------------------
% OUTPUT                                                                                                   
% -------------------------------------------------------------------------
% LDs: structure containing information about local datasets (LDs)
% and local models
%
%	LDs.X{i}:  matrix of regressors belonging to the i-th local dataset
%	(each row is a point).
%   LDs.y{i}: vector of outputs belonging to the i-th local dataset.
%   LDs.pos{i}: vector of indexes of datapoints in the i-th local
%   dataset, e.g. Xid(LDs.pos{i}(1),:) is the first regressor in the i-th
%   LD. 
%   LDs.weights{i} weight associated to the i-th datapoint used in
%   weighted LS for computing mode PVs. 
%   LDs.meanX{i} average of regressors in the i-th local dataset. 
%   LDs.models{i} parameters of the i-th local model.
%   LDs.models_var{i} INVERSE variance of the i-th local model.
%   LDs.X_ivar{i} INVERSE of the variance of the the regressors in the
%   i-th LD.

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

global idpar;
if ~isstruct(idpar),
    hit_error;
end

% Extract c
c=idpar.c; 

% ndata: # of data
% ndim: dimension of the X domain
ndata=size(Xid,1);
ndim=size(Xid,2);



% Build local datasets in the X space
% -------------------------------------------------------------

% LDs: structure containing all the informations about LDs and local models
%
%	  LDs.X: cell array: LDs.X{i}  matrix of X-points belonging to the i-th local dataset
%						 (each row is a point)
%     LDs.y: cell array: LDs.y{i} vector of y-points corresponding to the X-points
%						 belonging to the i-th local dataset.
%     LDs.pos: cell array: LDs.pos{i} vector of indices of the X points in the i-th local dataset.
%						 Such indices identify the position (i.e. the # of the row)
%						 in the matrix Xid of each point
%     LDs.weights: cell array: LDs.weights{i} =weight associated to y(i)
%						 Such weights should be used in the weighted LS.
%						 They can be produced for instance, when iterating
%						 the primal - extended dual identification
%     LDs.meanX: cell array: LDs.meanX{i} average of the X points in the i-th local dataset
%
%     LDs.models: cell array: LDs.models{i} parameters of the i-th local model
%
%     LDs.models_ivar: cell array: LDs.models_ivar{i} INVERSE of the variance of the i-th local model
%
%     LDs.X_ivar: cell array: LDs.X_var{i} INVERSE of the variance of the the X-points in the i-th LD
%
% Note that the number of local datasets equals the number of datapoints


% cand and Xdif are temporary vectors

cand=[];
Xdif=[];

LDs.X=cell(ndata,1);
LDs.y=cell(ndata,1);
LDs.pos=cell(ndata,1);
LDs.weights=cell(ndata,1);
LDs.meanX=cell(ndata,1);
LDs.models=cell(ndata,1);
LDs.models_ivar=cell(ndata,1);
LDs.X_ivar=cell(ndata,1);

for i=1:ndata
    cand=Xid(i,:);
    Xdif=Xid-ones(ndata,1)*cand;
    dist=[];
    for j=1:ndata
        dist=[dist;norm(Xdif(j,:),2)];  %compute the distances
    end

    % temp is temporary vector storing a set of indices
    temp=1:ndata;
    % increase dist with the indices in the first column
    dist=[temp',dist];
    % order the distances and the indices in one shot!
    dist=sortrows(dist,2);

    % tempXcl, tempycl tempindex and tempweight are temporary vectors
    % storing the informations about the i-the cluster

    % by convection, the first X-point of each LD
    % is the labeling point of the LD

    tempXcl=[Xid(i,:)];
    tempycl=[yid(i)];
    tempindex=[i];


    % assign the weight corresponding to the first point of the LD

    tempweight=Weight_primal(i);


    % add the remaining quantities for all the other points in the cluster
    % NB: the vector dist(:,1) stores the indices  in the X and y matrices
    % of the points belonging to the i-th cluster

    for j=2:c
        tempXcl=[tempXcl;Xid(dist(j,1),:)];
        tempycl=[tempycl;yid(dist(j,1))];
        tempindex=[tempindex,dist(j,1)];
        tempweight=[tempweight;Weight_primal(dist(j,1))];
    end

    % store the informations in the cell array LDs
    LDs.X{i}=tempXcl;
    LDs.y{i}=tempycl(:);
    LDs.pos{i}=tempindex(:);
    LDs.weights{i}=tempweight(:);
    LDs.meanX{i}=mean(tempXcl)'; % is based on the normalized data !!
end

% Identification of each linear model associated to a cluster
% -------------------------------------------------------------

% LDs.models:  stores the parameters of the i-th local model, i.e.
%		   	model identified by using the data of the i-th cluster.
%         LDs.models{i} is computed using WLS (or LS) on the normalized
%         datapoints


for i=1:ndata
    % add a vector of ones for estimating the constant through WLS
    Phi_temp=[LDs.X{i},ones(c,1)];
    [temp,stdx,mse,S]=hit_lscov(Phi_temp,LDs.y{i},LDs.weights{i});
    LDs.models{i}=temp;
    LDs.models_ivar{i}=inv(S);
    LDs.X_ivar{i}=inv(cov(LDs.X{i})); 
    if ~all(all(isfinite(LDs.models_ivar{i}))) | ~all(all(isfinite(LDs.X_ivar{i})))
        fprintf('\n ======================================================== ');
        fprintf('\n ERROR IN CREATING LDs: some covariance matrices of local'); 
        fprintf('\n models do not have finite entries.');
        fprintf('\n A possible cause is that data are almost noiseless and');
        fprintf('\n local models are perfectly linear.');
        fprintf('\n Rough remedies: ');
        fprintf('\n 1) increase idpar.c;');
        fprintf('\n 2) add a little noise to output data.');
        fprintf('\n hit_regression.m ends here. ');
        fprintf('\n ======================================================== \n');
        error('hit_regression:end','hit_regression ends here.')
    end
end






