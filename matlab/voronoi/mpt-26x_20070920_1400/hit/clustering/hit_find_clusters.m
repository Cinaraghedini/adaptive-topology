function [new_centers,new_Inl,new_class,err,clustvalid]=hit_find_clusters(FV,s,dim_data)
%HIT_FIND_CLUSTERS Interface with clustering algorithms.
%
% -------------------------------------------------------------------------
% DESCRIPTION
% -------------------------------------------------------------------------
% [new_centers,new_Inl,new_class,err,clustvalid]= ...
% hit_find_clusters(FV,s,dim_data)
%
% -------------------------------------------------------------------------
% INPUT
% -------------------------------------------------------------------------
% FV(i,:): i-th xi-point to be clustered.
% s: number of clusters to be found (used just by algorithms that do not
% perform the estimation of the cluster number).
% dim_data: dimension of the xi-points.
%
% -------------------------------------------------------------------------
% OUTPUT                                                                                                   
% -------------------------------------------------------------------------
% new_centers(i,:): center of the i-th cluster (it is a ROW vector).
% new_Inl: indexes of datapoints that are inliers after clustering (i.e.
% not discarded by the clustering algorithm).
% new_class(i): classification of the i-th inlier. new_class(i)=j means
% that the i-th inlier elongs to the j-th cluster.
% err: value of the clustering cost functional at the optimal solution.
% clustvalid: structure containing information about the
% clustering.results. It is equal to the structure 'spec' produced by
% specific clustering algorithms such as hit_kmeans or hit_single_linkage.

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

% Initialization: Inliers=all points

Inl=[FV.points{:}]';
w_Inl=[FV.weights{:}]';
IR_Inl=FV.IR;
nInl=size(Inl,1);


switch upper(idpar.clustalgo.name)
case{'KMEANS'} 
    % Runs the modified K-means algorithm.
    % Note that the algorithm is called even if s=1 (one cluster)
    % This allows to set correctly the outputs of hit_find_clusters
    % and to have useful cluster validation measures
    % However, silhouette values are computed only if s>1
    
    % Set the options vector for K-means (see kmeans_mod2 for details) 
    options(1)  = 0;		% Do not print out error values.
    options(2) = 1e-4;
    options(3)=1e-4;
    options(4) = 100;   % Max n. of cycles in k-means
    
    % Set the options of the modified K-means
    % according to what the used specified
    
    % fields checking 
    
    opt.init_centers=idpar.clustalgo.kmeans.init_centers;
    opt.centers=idpar.clustalgo.kmeans.centers;
    opt.repetitions=idpar.clustalgo.kmeans.repetitions;
    opt.s=s;
    switch lower(idpar.clustalgo.kmeans.init_centers)
    case{'covariances'}
        opt.IR=IR_Inl;
    case{'scalars'}
        opt.w=w_Inl;
    end
    switch lower(idpar.clustalgo.kmeans.centers)
    case{'covariances'}
        opt.IR=IR_Inl;
        case{'scalars'}
        opt.w=w_Inl;
        
    end
    

[new_centers,err,new_Inl,new_class,etime,kmeans_stat] = hit_kmeans(Inl,opt);
    clustvalid=kmeans_stat;
case{'SL'}
    % Compute the clusters by using the single-linkage algorithm 
    %		opt: structure of specific parameters
	%            MANDATORY FIELDS
	%            opt.guess: guessed minimal distanceS
    %            opt.plot_step: 'Y', 'N' )
        
    opt.guess=idpar.clustalgo.sl.guess_min_dist;
    opt.plot_steps=idpar.clustalgo.sl.plot_steps;
    opt.plot_fig=idpar.clustalgo.sl.plot_fig;
    % add one cluster for the outliers !
    opt.s=s+1;
    [new_centers,err,new_Inl,new_class,etime,sl_stats] = hit_single_linkage(Inl, opt);
    clustvalid=sl_stats;
	%             spec.costs: spec.costs(i) is the costs of the i-th cluster
	%             spec.n_mclust: number of  clusters
    % note that all the points in Inl belong to some micro-cluster
end
