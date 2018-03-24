function [topreserve,isequalto]=hit_remove_duplicates(LDs,xi)
%HIT_REMOVE_DUPLICATES Remove duplicates of xi-points before clustering.
%
% -------------------------------------------------------------------------
% DESCRIPTION
% -------------------------------------------------------------------------
%
% [topreserve,isequalto]=hit_remove_duplicates(LDs,xi)
% 
% -------------------------------------------------------------------------
% INPUT
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
%   i-th LD.Xid: matrix containing the datapoints in the regressor set.
%       Each row is a datapoint.
%
% xi: structure containing information about the xi-points (they are either
% FVs or LPVs).
%
%	xi.points{i}: xi-point based on the i-th LDs.
%   xi.IR{i}: INVERSE of the covariance of the i-th xi-point.
%   xi.weights{i}: scalar confidence measure of the i-th feature vector.
%
% -------------------------------------------------------------------------
% OUTPUT                                                                                                   
% -------------------------------------------------------------------------
% topreserve: indexes of xi-points that must be preserved.
% isequalto{i}: indexes of xi-points that are equal to the point indexed by
% topreserve(i). It can be an empty vector (if a point does not have
% duplicates)

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

nLDs=length(LDs.y);
isequalto=cell(nLDs,1);

toremove=[];
for i=1:nLDs
    for j=i+1:nLDs
        if isempty(setxor(LDs.pos{i},LDs.pos{j}))
            toremove=[toremove,i];
            %isequalto{i}=[isequalto{i};j];
            isequalto{j}=[isequalto{j};isequalto{i};i];
            isequalto{j}=unique(isequalto{j});
            %delete useless dependecies
            isequalto{i}=[];
        end
    end
end
topreserve=setdiff([1:nLDs],toremove);
isequalto=isequalto(topreserve);