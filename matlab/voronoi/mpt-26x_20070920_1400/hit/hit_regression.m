function [idmodes,F,xi,LDs,inliers]=hit_regression(Xid,yid)
%HIT_REGRESSION Main function for PWA regression.
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
% -------------------------------------------------------------------------
% OUTPUT                                                                                                   
% -------------------------------------------------------------------------
% idmodes: structure containing all information on the identified models.
%
%   idmodes.par{i}: is the PVs of the i-th mode.
%   idmodes.regions(i): is the region of the i-th mode. Regions are
%   subsets of idpar.Regressor_set. 
%   idmodes.cov{i}; is the covariance of the PV of the i-th mode. 
%   idmodes.Regressor_set: automatically computed polytope defining the
%   regressor set. This field is present only if idpar.Regressor_set is
%   empty or nonexistent.    
%   idmodes.s: meaningful number of modes estimated during regression (it
%   might differ from idpar.s).
%   idmodes.regions_sim(i): is the region for simulating the i-th mode.
%   These regions are subsets of idpar.Regressor_set_sim and are defined
%   only if idpar.Regressor_set_sim is defined.
%   idmodes.stat_reattr: fraction of points reattribted to modes, if
%   reattribution has been performed.
%   idmodes.clust_valid: structure containing information about the
%   clustering results.
%   idmodes.pattern_rec_valid: structure containing information about the
%   pattern recognition results.
%   idmodes.Weight_primal: weights used in creating LDs (usually it is a
%   vector of 1's) 
%   idmodes.adjacences: each row is a pair (i,j), i<j indicating that the
%   regions i and j are adjacent The constraint i<j avoids storing both
%   pairs (i,j) and (j,i).
% F: structure containing information about the mode datasets (the
% classified datapoints that are also inliers, i.e. not discarded during
% regression).
%
%   F.X{i}: matrix of regressors assigned to the i-th mode dataset.
%   F.y: cell array: F.y{i} vector of outputs assigned to the i-th mode
%   dataset. 
%   F.pos{i}: indexes of the points assigned to the i-th mode dataset. 
%   F.pos{i}(1)=5 means that Xid(5,:) and yid(5) are the first points
%   composing the i-th mode dataset.
%
% xi: structure containing information about the xi-points (they are either
% FVs or LPVs).
%
%	xi.points{i}: xi-point based on the i-th LDs.
%   xi.IR{i}: INVERSE of the covariance of the i-th xi-point.
%   xi.weights{i}: scalar confidence measure of the i-th feature vector.
%
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
%
% inliers: structure containing information on the inliers
%
%	inliers.pos(j): index of the j-th inlier in Xid and yid
%	(i.e.Xid(inliers.pos(j),:) and yid(inliers.pos(j)) are the j-th
%	inliers).
%   inliers.class(j): classification of the j-th inlier.

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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHECK DATA AND FIELDS OF IDPAR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check if data are well-defined

if ~all(all(isfinite(Xid)))|~all(isfinite(yid))
    fprintf('\n ======================================================== ');
    fprintf('\n DATA CHECK ERROR: some regressors or outputs.');
    fprintf('\n contains non-numeric entries (e.g. inf or nan). ');
    fprintf('\n ======================================================== \n');
    error('hit_regression:end','hit_regression ends here.')
end

% ndata: # of datapoints
% ndim: dimension of the X-space
ndim=size(Xid,2);
ndata=size(Xid,1);

if ndata<=ndim
    fprintf('\n ======================================================== ');
    fprintf('\n DATA CHECK ERROR: insufficient number of datapoints.');
    fprintf('\n Check if the regressors are stored as ROWS of Xid. ');
    fprintf('\n ======================================================== \n');
    error('hit_regression:end','hit_regression ends here.')
end

% put the outputs in a column vector

yid=yid(:);

% checks mandatory fields for building LDs

if ~isfield(idpar,'c') | isempty(idpar.c)
    fprintf('\n ======================================================== ');
    fprintf('\n IDPAR CHECK ERROR: no size of LDs specified.');
    fprintf('\n Set idpar.c and retry.');
    fprintf('\n ======================================================== \n');
    error('hit_regression:end','hit_regression ends here.')
elseif idpar.c<2 
        % this check i needed because hit_create_LDs assumes that idpar.c>=2
        % Moreover, idpar.c<2 has little meaning even if the regressors are
        % scalars
        fprintf('\n ======================================================== ');
        fprintf('\n IDPAR CHECK ERROR: idpar.c is less than 2.');
        fprintf('\n Increase idpar.c and retry.');
        fprintf('\n ======================================================== \n');
        error('hit_regression:end','hit_regression ends here.')
end

% Check if there are enough data to have at least one LD
if ndata<=idpar.c
    fprintf('\n ======================================================== ');
    fprintf('\n DATA CHECK ERROR: insufficient number of datapoints.');
    fprintf('\n The number of datapoints is less than idpar.c and no');
    fprintf('\n local dataset can be created.');
    fprintf('\n Hint: check if the regressors are stored as ROWS of Xid.');
    fprintf('\n ======================================================== \n');
    error('hit_regression:end','hit_regression ends here.')
end
% checks mandatory fields for clustering

if ~isfield(idpar.clustalgo,'name') | ~iscellstr({idpar.clustalgo.name})
    fprintf('\n ======================================================== ');
    fprintf('\n IDPAR CHECK ERROR: no clustering algorithm specified');
    fprintf('\n Set idpar.clustalgo.name and retry.');
    fprintf('\n ======================================================== \n');
    error('hit_regression:end','hit_regression ends here.')
else
    % check the field values
    switch lower(idpar.clustalgo.name)
        case{'kmeans','sl'}
        otherwise
            fprintf('\n ======================================================== ');
            fprintf('\n IDPAR CHECK ERROR: the clustering algorithm specified');
            fprintf('\n does not match one of those available in the HIT toolbox ');
            fprintf('\n check idpar.clustalgo.name.');
            fprintf('\n ======================================================== \n');
            error('hit_regression:end','hit_regression ends here.')
    end
end
if ~isfield(idpar,'what_to_clust') | ~iscellstr({idpar.what_to_clust})
    fprintf('\n ======================================================== ');
    fprintf('\n IDPAR CHECK ERROR: cluster LPVs or FVs ?');
    fprintf('\n Set idpar.what_to_clust and retry.');
    fprintf('\n ======================================================== \n');
    error('hit_regression:end','hit_regression ends here.')
else
    % check the field values
    switch lower(idpar.what_to_clust)
        case{'lpvs','fvs'}
        otherwise
            fprintf('\n ======================================================== ');
            fprintf('\n IDPAR CHECK ERROR: the value of idpar.what_to_clust');
            fprintf('\n does not match one of those allowed in the HIT toolbox.');
            fprintf('\n ======================================================== \n');
            error('hit_regression:end','hit_regression ends here.')
    end
end
if ~isfield(idpar.clustalgo,'remove_duplicates') | ~iscellstr({idpar.clustalgo.remove_duplicates})
    fprintf('\n ======================================================== ');
    fprintf('\n IDPAR CHECK ERROR: remove duplicates before clustering ?');
    fprintf('\n Set idpar.clustalgo.remove_duplicate and retry.');
    fprintf('\n ======================================================== \n');
    error('hit_regression:end','hit_regression ends here.')
else
    % check the field values
    switch lower(idpar.clustalgo.remove_duplicates)
        case{'y','n'}
        otherwise
            fprintf('\n ======================================================== ');
            fprintf('\n IDPAR CHECK ERROR: the value of');
            fprintf('\n idpar.clustalgo.remove_duplicates must be Y or N.');
            fprintf('\n ======================================================== \n');
            error('hit_regression:end','hit_regression ends here.')
    end
end


% check the parameters needed by each clustering algorithm.
% Some algorithms require also idpar.s as input
% For instance, kmeans need idpar.s, SL does not need it

switch lower(idpar.clustalgo.name)
    case{'kmeans'}
        if ~isfield(idpar,'s') | isempty(idpar.s)
            fprintf('\n ======================================================== ');
            fprintf('\n IDPAR CHECK ERROR: no mode number specified but the mode');
            fprintf('\n number is needed by the clustering algorithm chosen.');
            fprintf('\n Set idpar.s and retry');
            fprintf('\n ======================================================== \n');
            error('hit_regression:end','hit_regression ends here.')
        end
        if min(idpar.s)<1
            fprintf('\n ======================================================== ');
            fprintf('\n IDPAR CHECK ERROR: a mode number <1 has been set');
            fprintf('\n but the entries of idpar.s must be >= 1.');
            fprintf('\n Set idpar.s and retry');
            fprintf('\n ======================================================== \n');
            error('hit_regression:end','hit_regression ends here.')
        end
        % make sure that idpar.s is a row
        idpar.s=idpar.s(:)';
        if ~isfield(idpar.clustalgo.kmeans,'init_centers') | ~iscellstr({idpar.clustalgo.kmeans.init_centers})
            fprintf('\n ======================================================== ');
            fprintf('\n IDPAR CHECK ERROR: no algorithm has been chosen for');
            fprintf('\n initilizing the cluster centers in kmeans.');
            fprintf('\n Check idpar.clustalgo.kmeans.init_centers.');
            fprintf('\n ======================================================== \n');
            error('hit_regression:end','hit_regression ends here.')
        else
            switch lower(idpar.clustalgo.kmeans.init_centers)
                case{'covariances','scalars'}
                otherwise
                    fprintf('\n ======================================================== ');
                    fprintf('\n IDPAR CHECK ERROR: no valid algorithm');
                    fprintf('\n for initializing the centers in kmeans');
                    fprintf('\n has been specified.');
                    fprintf('\n Check idpar.clustalgo.kmeans.init_centers.');
                    fprintf('\n ======================================================== \n');
                    error('hit_regression:end','hit_regression ends here.')
            end
        end
        if ~isfield(idpar.clustalgo.kmeans,'centers') | ~iscellstr({idpar.clustalgo.kmeans.centers})
            fprintf('\n ======================================================== ');
            fprintf('\n IDPAR CHECK ERROR: no algorithm has been chosen for');
            fprintf('\n computing the cluster centers in kmeans.');
            fprintf('\n Check idpar.clustalgo.kmeans.centers.');
            fprintf('\n ======================================================== \n');
            error('hit_regression:end','hit_regression ends here.')
        else
            switch lower(idpar.clustalgo.kmeans.centers)
                case{'covariances','scalars'}
                otherwise
                    fprintf('\n ======================================================== ');
                    fprintf('\n IDPAR CHECK ERROR: no valid algorithm');
                    fprintf('\n for updating the centers in kmeans');
                    fprintf('\n has been specified.');
                    fprintf('\n Check idpar.clustalgo.kmeans.centers.');
                    fprintf('\n ======================================================== \n');
                    error('hit_regression:end','hit_regression ends here.')
            end
        end
        if ~isfield(idpar.clustalgo.kmeans,'repetitions') | isempty(idpar.clustalgo.kmeans.repetitions)
            fprintf('\n ======================================================== ');
            fprintf('\n IDPAR CHECK ERROR: the number of times');
            fprintf('\n kmeans must be run with different initializations');
            fprintf('\n has not been specified.');
            fprintf('\n Check idpar.clustalgo.kmeans.repetitions.');
            fprintf('\n ======================================================== \n');
            error('hit_regression:end','hit_regression ends here.')
        end
    case{'sl'}
        if ~isfield(idpar.clustalgo.sl,'guess_min_dist') | isempty(idpar.clustalgo.sl.guess_min_dist)
            fprintf('\n ======================================================== ');
            fprintf('\n IDPAR CHECK ERROR: no minimal distance between');
            fprintf('\n clusters has been specified for single-linkage clustering.');
            fprintf('\n Check idpar.clustalgo.SL.guess_min_dist.');
            fprintf('\n ======================================================== \n');
            error('hit_regression:end','hit_regression ends here.')
        end
        if ~isfield(idpar.clustalgo.sl,'plot_steps') | ~iscellstr({idpar.clustalgo.sl.plot_steps})
            fprintf('\n ======================================================== ');
            fprintf('\n IDPAR CHECK ERROR: plot each step of single-linkage ?');
            fprintf('\n Check idpar.clustalgo.SL.plot_steps.');
            fprintf('\n ======================================================== \n');
            error('hit_regression:end','hit_regression ends here.')
        else
            % check the field values
            switch lower(idpar.clustalgo.sl.plot_steps)
                case{'y','n'}
                otherwise
                    fprintf('\n ======================================================== ');
                    fprintf('\n IDPAR CHECK ERROR: the value of');
                    fprintf('\n idpar.clustalgo.sl.plot_step must be Y or N.');
                    fprintf('\n ======================================================== \n');
                    error('hit_regression:end','hit_regression ends here.')
            end
        end
        if strcmp(lower(idpar.clustalgo.sl.plot_steps),'y')
            if ~isfield(idpar.clustalgo.sl,'plot_fig') | isempty(idpar.clustalgo.sl.plot_fig)
                fprintf('\n ======================================================== ');
                fprintf('\n IDPAR CHECK ERROR: no figure number has been specified');
                fprintf('\n for plotting steps of single-linkage clustering.');
                fprintf('\n Check idpar.clustalgo.sl.plot_fig.');
                fprintf('\n ======================================================== \n');
                error('hit_regression:end','hit_regression ends here.')
            end
        end
end



% Check parameters needed in the clustering post-processing
% and for finding the modes
if ~isfield(idpar,'mix_detect') | ~iscellstr({idpar.mix_detect})
    fprintf('\n ======================================================== ');
    fprintf('\n IDPAR CHECK ERROR: detect a posteriori possible mixed datapoints ?');
    fprintf('\n Set idpar.mix_detect and retry.');
    fprintf('\n ======================================================== \n');
    error('hit_regression:end','hit_regression ends here.')
else
    % check the field values
    switch lower(idpar.mix_detect)
        case{'y','n'}
        otherwise
            fprintf('\n ======================================================== ');
            fprintf('\n IDPAR CHECK ERROR: the value of');
            fprintf('\n idpar.clustalgo.remove_duplicates must be Y or N.');
            fprintf('\n ======================================================== \n');
            error('hit_regression:end','hit_regression ends here.')
    end
end
if ~isfield(idpar,'continuity') | ~iscellstr({idpar.continuity})
    fprintf('\n ======================================================== ');
    fprintf('\n IDPAR CHECK ERROR: the model is continuous or not ?');
    fprintf('\n Set idpar.continuity to C or D.');
    fprintf('\n ======================================================== \n');
    error('hit_regression:end','hit_regression ends here.')
else
    % check the field values
    switch lower(idpar.continuity)
        case{'c','d'}
        otherwise
            fprintf('\n ======================================================== ');
            fprintf('\n IDPAR CHECK ERROR: the value of');
            fprintf('\n idpar.continuity must be C or D.');
            fprintf('\n ======================================================== \n');
            error('hit_regression:end','hit_regression ends here.')
    end
end
if ~isfield(idpar,'discard_threshold_factor') | isempty(idpar.discard_threshold_factor)
    fprintf('\n ======================================================== ');
    fprintf('\n IDPAR CHECK ERROR: the discard factor for removing modes');
    fprintf('\n with few points has not been defined');
    fprintf('\n Set idpar.discard_threshold_factor and retry.');
    fprintf('\n ======================================================== \n');
    error('hit_regression:end','hit_regression ends here.')
end
if ~isfield(idpar,'patt_rec_algo') | ~iscellstr({idpar.patt_rec_algo})
    fprintf('\n ======================================================== ');
    fprintf('\n IDPAR CHECK ERROR: the pattern recognition algorithm');
    fprintf('\n (for finding the mode regions) is not specified.');
    fprintf('\n Set idpar.patt_rec_algo.');
    fprintf('\n ======================================================== \n');
    error('hit_regression:end','hit_regression ends here.')
else
    % check the field values
    switch lower(idpar.patt_rec_algo)
        case{'svc','psvc','mrlp'}
        otherwise
            fprintf('\n ======================================================== ');
            fprintf('\n IDPAR CHECK ERROR: the value of idpar.patt_rec_algo');
            fprintf('\n does not match one of the pattern recognition methods');
            fprintf('\n available in the HIT toolbox.');
            fprintf('\n ======================================================== \n');
            error('hit_regression:end','hit_regression ends here.')
    end
end
% Fields to be set automatically if not specified by the user
% They are saved into idmodes to let the user keep track of
% the automatic choices

% Set automatically the regressor set if it has not been supplied by the
% user, and update idmodes accordingly
if ~isfield(idpar,'Regressor_set') | isempty(idpar.Regressor_set)
    Regressor_set=polytope([eye(ndim);-eye(ndim)],[max(Xid)';-min(Xid)']);
    idmodes.Regressor_set=Regressor_set;
else
    [BigA,BigB]=double(idpar.Regressor_set);
    % Check if the dimension of idpar.Regressor_set matches the dimension of the
    % regressors
    if size(BigA,2)~=ndim
        fprintf('\n ======================================================== ');
        fprintf('\n DATA CHECK ERROR: the dimension of idpar.Regressor_set');
        fprintf('\n does not match the dimension of regressors. ');
        fprintf('\n ======================================================== \n');
        error('hit_regression:end','hit_regression ends here.')
    end
    Regressor_set=idpar.Regressor_set;
end
% Check if the dimension of idpar.Regressor_set_sim matches the dimension of the
% regressors
if isfield(idpar,'Regressor_set_sim') & ~isempty(idpar.Regressor_set_sim)
    [BigA,BigB]=double(idpar.Regressor_set_sim);
    if size(BigA,2)~=ndim
        fprintf('\n ======================================================== ');
        fprintf('\n DATA CHECK ERROR: the dimension of ');
        fprintf('\n idpar.Regressor_set_sim does not match the dimension ');
        fprintf('\n of regressors. ');
        fprintf('\n ======================================================== \n');
        error('hit_regression:end','hit_regression ends here.')
    end
end
% Set automatically the the a priori weights on data, if they have
% been not subblied by the user, and update idmodes accordingly
if ~isfield(idpar,'Weight_primal') | isempty(idpar.Weight_primal)
    Weight_primal=ones(size(yid)); % no weights associated to data points
    % Store it for the user
    idmodes.Weight_primal=Weight_primal;
else
    Weight_primal=idpar.Weight_primal;
end

% Fields to be set automatically (if not specified by the user)
% and saved into idpar to let the user keep track of
% the automatic choices

% Set automatically to 'n' the questions if no indication
% has been given by the user and update idpar accordingly
% this is the only field of idpar that can be changed by hit_regression

if ~isfield(idpar,'YNquestions')
    YNquestions='n';
else
    YNquestions=idpar.YNquestions;
end

% check the fields needed in plotpar for the correct running of 
% hit_regression and the functions called 
hit_plotpar_check({'marker_modes','plot_xi_points_yn','plot_xi_points_fig','plot_class_xi_points_yn','plot_class_xi_points_fig'})

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END OF CHECKING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



fprintf('\n *****************************');
fprintf('\n *** PWA regression begins ...');
fprintf('\n ****************************** \n');

% CREATE LOCAL models (LDs)
% -------------------------------------------------------------

%  LDs: structure containing all the information about LDs and local models, see 
%  hit_create_LDs.m
%

fprintf('\n --- CREATING LOCAL MODELS ...');

LDs=hit_create_LDs(Xid,yid,Weight_primal);


%warning off;




% BUILD THE xi VECTORS
% these are LPVs or FVs according to idpar.what_to_clust
% -------------------------------------------------------------


% xi: structure containing all the informations about the xi vectors
%
%	  xi.points: cell array: xi.points{i} si the xi-vector based on the i-th LD
%     xi.IR: cell array: xi.IR{i} INVERSE of the covariance of the i-th xi-vector
%     xi.weights: cell array: xi.weights{i} store the scalar confidence measure of the i-th xi-vector

xi.points=cell(ndata,1);
xi.IR=cell(ndata,1);
xi.weights=cell(ndata,1);

for i=1:ndata
    switch lower(idpar.what_to_clust)
        case{'fvs'}
            xi.points{i}=[LDs.models{i};LDs.meanX{i}]; %column vector
            xi.IR{i}=[LDs.models_ivar{i},zeros(ndim+1,ndim);zeros(ndim,ndim+1), LDs.X_ivar{i}];
        case{'lpvs'}
            xi.points{i}=[LDs.models{i}]; %column vector
            xi.IR{i}=LDs.models_ivar{i};
    end
end

% dimxi: dimension of each vector in the feature space
dimxi = size(xi.points{1},1);

for i=1:ndata
    xi.weights{i}=1/sqrt((2*pi)^dimxi*inv(det(xi.IR{i})));
end
if ~all(isfinite([xi.weights{:}]))
        fprintf('\n ======================================================== ');
        fprintf('\n ERROR IN COMPUTING THE SCALAR WEIGHTS TO xi-POINTS:'); 
        fprintf('\n some weights are not finite.');
        fprintf('\n A possible cause is that data are almost noiseless.');
        fprintf('\n Rough remedy: add a little noise to output data.');
        fprintf('\n hit_regression.m ends here. ');
        fprintf('\n ======================================================== \n');
        error('hit_regression:end','hit_regression ends here.')
end


% Place here the automatic computation/suggestion of some
% parameters specific to the clustering algorithm chosen.
% At present, no parameter is suggested.
% The next commented lines give an example when choosing the
% algorithm SL.
%
% switch lower(idpar.clustalgo.name)
%     case{'SL'}
%         switch lower(idpar.clustalgo.suggest_min_dist)
%             case{'y'} %suggest a minimal distance
%                 % PLACE YOUR CODE HERE FOR SUGGESTING A MINIMAL DISTANCE
%                 % BETWEEN CLUSTERS
%         end
% end

% PLOT THE FIRST THREE COMPONENTS OF THE FVs/LPVs
% BEFORE CLUSTERING

% This figure is meaningful just if the FVs/LPVs are of dimension
% 2 or 3. In this case, clusters should be
% seen by visual inspection. If this does not happens,
% it is recommended to lower the parameter idpar.c.
% However, if the dimension is >3, just the projections of xi on
% the first three dimensions are plotted. Hence, clusters
% may be not visible, even if they are well-separated.
% If clusters can be seen visually, their number
% corresponds to the number of modes. This can be used to set
% the parameter idpar.s
if strcmp(lower(plotpar.plot_xi_points_yn),'y')
hit_plot_clust([xi.points{:}]',[],[1:ndata]',ones(ndata,1),['xb'],plotpar.plot_xi_points_fig);
switch lower(idpar.what_to_clust)
    case{'fvs'}
        title('First three componenents of the FVs')
    case{'lpvs'}
        title('First 2 (or 3) componenents of the LPVs')
end
hold off
end


%
% CLUSTERING
% at the end, the number of modes will be stored in the variable s
%

fprintf('\n --- CLUSTERING xi-points ...');
% set the mode number s, if needed, otherwise
        % set s=0
        switch lower(idpar.clustalgo.name)
            case{'kmeans'}
                s=idpar.s;
            case{'sl'}
                % s can be whatever value: SL finds automatically s
                s=0;
        end


% remove duplicates among xi points if required according to
% idpar.clustalgo.remove_duplicates
switch upper(idpar.clustalgo.remove_duplicates)
    case{'Y'}
        fprintf('\n Removing duplicate xi-points before clustering ...')
        [topreserve,isequalto]=hit_remove_duplicates(LDs,xi);
        xi1.points={xi.points{topreserve}};
        xi1.IR={xi.IR{topreserve}};
        xi1.weights={xi.weights{topreserve}};

        % Call the clustering algorithm
        [centers,inliers.pos,inliers.class,clust_ssr]=hit_find_clusters(xi1,s,ndim);
        % associate the duplicates to the right clusters
        % and store the results in ind_inl and ind_class_inl
        ind_inl_new=[];
        ind_class_inl_new=[];
        for kkk=1:length(ind_inl)
            ind_inl_new=[ind_inl_new;topreserve(inliers.pos(kkk))];
            ind_class_inl_new=[ind_class_inl_new; inliers.class(kkk)];
            sz=length(isequalto{inliers.pos(kkk)});
            if sz>0;
                ind_inl_new=[ind_inl_new;isequalto{inliers.pos(kkk)}];
                ind_class_inl_new=[ind_class_inl_new;inliers.class(kkk)*ones(sz,1)];
            end
        end
        inliers.pos=ind_inl_new;
        inliers.class=ind_class_inl_new;
    case{'N'}
        % Call the clustering algorithm
        [centers,inliers.pos,inliers.class,clust_ssr,clust_valid]=hit_find_clusters(xi,s,ndim);
end %ends switch upper(idpar.clustalgo.remove_duplicate_FVs)

% update (or define) the number of modes
s=size(centers,1);


% save the clustering validation indicators in idmodes  
idmodes.clust_valid=clust_valid;
% Plot the clusters of xi-points
% the routine plot_clust3d extract, from vectors and centers
% in the feature space, the first three components of the
% feature vectors (2 components if the FVs are 2-dimensional)
%
% The plot is meaningful for checking if the estimated centers
% describe the clusters in the feature space
% -------------------------------------------------------------
if strcmp(lower(plotpar.plot_class_xi_points_yn),'y')
hit_plot_clust([xi.points{:}]',centers,inliers.pos,inliers.class,plotpar.marker_modes,plotpar.plot_class_xi_points_fig)
switch lower(idpar.what_to_clust)
    case{'fvs'}
        title('First three componenents of the classified FVs')
    case{'lpvs'}
        title('First 2 (or 3) componenents of the classified LPVs')
end
end

% CLUSTERING - POST PROCESSING


% Update s if some clustering algorithm
% estimating the number of models has been used

switch lower(idpar.clustalgo.name)
    case{'sl'}
        s=size(centers,1);
        fprintf('\n Number of modes estimated by Single-Linkage clustering: %d',s);
end

% Remove clusters that have less than 2*ndim+1 data points
% This points will be removed from the inliers and the number of modes
% is reduced accordingly

% Remove clusters that have less than idpar.c points and decrement
% the number of modes accordingly.
% The points in such clusters are removed from the list of inliers

% These tests have also the effect of removing empty clusters, 
% a phenomenon that can happen in clustering (e.g. in kmeans)
[inliers,s]=hit_check_inliers(inliers,s,ndim);

switch lower(idpar.YNquestions)
    case{'y'}
        contYN=input('Continue and find the submodels (Y/N) ?','s')
    case{'n'}
        contYN='y';
end
switch lower(contYN)
    case{'n'}
        return
end

%
% MODE ESTIMATION
%
        %%%% TO DO
        % hit_reattribute_suspected
        % is a batch files ! It should be transformed into a function.
        % It just modifies inliers and creates exchange and stat_reattr
        %%%% END TO DO
        %%%% TO DO
        % hit_attribute_outliers
        % is a batch files ! It should be transformed into a function.
        % It just modify inliers and set the flag outl_attributed=1, if
        % some reattribution took place
        %%%% END TO DO
        
% list of outlier indexes
out_pos=setdiff([1:ndata],inliers.pos);

if isempty(out_pos) %there is no outlier
    switch upper(idpar.mix_detect)
        case{'Y'} % reattribute xi-points
            % compute mode PVs and regions of a discontinuous PWA model
            % NB: also the regions are needed for reattributing xi-points
            continuity='d';
            fprintf('\n --- COMPUTING MODES ...');
            [F,idmodes1,adjacences]=...
                hit_compute_modes(Xid,yid,ndata,s,xi,inliers,Regressor_set,continuity);
            % reattribute xi-points
            fprintf('\n --- REATTRIBUTING DATAPOINTS SUSPECTED TO BE MIXED ...');
            stat_reattr=0; % default statistic of reattribution = 0
            hit_reattribute_suspected
            % Add the field stat_reattr to idmodes
            idmodes.stat_reattr=stat_reattr;
            if ~isempty(exchange) % some reattribution happened
                % update inliers and s
                [inliers,s]=hit_check_inliers(inliers,s,ndim);
            end
            % Compute the final modes only if some xi-points have been
            % reattributed or if the final PWA model is continuous
            % Otherwise, the final modes have been already computed few
            % lines above
            if ~isempty(exchange)| strcmp(idpar.continuity,'c')
                fprintf('\n --- COMPUTING MODES ...');
                [F,idmodes1,adjacences]=...
                    hit_compute_modes(Xid,yid,ndata,s,xi,inliers,Regressor_set,idpar.continuity);
            end
        case{'N'} % do not reattribute xi-points
            % Compute the final modes
            fprintf('\n --- COMPUTING MODES ...');
            [F,idmodes1,adjacences]=...
                hit_compute_modes(Xid,yid,ndata,s,xi,inliers,Regressor_set,idpar.continuity);
    end
else %there is some outlier
    switch upper(idpar.mix_detect)
        case{'Y'} % reattribute xi-points
            % compute mode PVs and regions of a discontinuous PWA model
            % NB: also the regions are needed for reattributing xi-points
            continuity='d';
            fprintf('\n --- COMPUTING MODES ...');
            [F,idmodes1,adjacences]=...
                hit_compute_modes(Xid,yid,ndata,s,xi,inliers,Regressor_set,continuity);
            % reattribute xi-points
            fprintf('\n --- REATTRIBUTING xi POINTS SUSPECTED TO BE MIXED ...');
            stat_reattr=0; % default statistic of reattribution = 0
            hit_reattribute_suspected
            % Add the field stat_reattr to idmodes
            idmodes.stat_reattr=stat_reattr;
            if ~isempty(exchange) % some reattribution happened
                % update inliers and s
                [inliers,s]=hit_check_inliers(inliers,s,ndim);
                %
                % MODE PVs ESTIMATION (discontinuous model)
                %
                % Here we need to update only idmodes1.par and not the regions
                % This is indeed the relevant information for the current
                % method used for attributing outliers
                % Moreover, the model is assumed to be discontinuous
                F=hit_compute_F(s,inliers,Xid,yid);
                fprintf('\n --- COMPUTING MODE PVs ...');
                [final_PVs,idmodes1.cov,...
                    Weight_primal]=hit_compute_disc_modes(F,ndata,[xi.weights{:}]);
                % stores the coefficients of the models in the structure
                % idmodes1
                for i=1:s
                    idmodes1.par{i}=(final_PVs(i,:))';
                end
            end
            % attribute outliers
            fprintf('\n --- ATTRIBUTING OUTLIERS TO MODES (MODES WILL BE RECOMPUTED) ...');
            hit_attribute_outliers
            % update inliers and s
            [inliers,s]=hit_check_inliers(inliers,s,ndim);
            % Compute the final modes because some outliers
            % have been reattributed for sure
            fprintf('\n --- COMPUTING MODES ...');
            [F,idmodes1,adjacences]=...
                hit_compute_modes(Xid,yid,ndata,s,xi,inliers,Regressor_set,idpar.continuity);
        case{'N'} % do not reattribute xi-points
            %
            % MODE PVs ESTIMATION (discontinuous model)
            %
            % Here we need to update only idmodes1.par and not the regions
            % This is indeed the relevant information for the current
            % method used for attributing outliers
            % Moreover, the model is aasumed to be discontinuous
            F=hit_compute_F(s,inliers,Xid,yid);
            fprintf('\n --- COMPUTING MODE PVs ...');
            [final_PVs,idmodes1.cov,...
                Weight_primal]=hit_compute_disc_modes(F,ndata,[xi.weights{:}]);
            % stores the coefficients of the models in the structure
            % idmodes1
            for i=1:s
                idmodes1.par{i}=(final_PVs(i,:))';
            end
            % attribute outliers
            fprintf('\n --- ATTRIBUTING OUTLIERS TO MODES (MODES WILL BE RECOMPUTED) ...');
            hit_attribute_outliers
            % update inliers and s
            [inliers,s]=hit_check_inliers(inliers,s,ndim);
            % COMPUTE FINAL MODES
            fprintf('\n --- COMPUTING MODES ...');
            [F,idmodes1,adjacences]=...
                hit_compute_modes(Xid,yid,ndata,s,xi,inliers,Regressor_set,idpar.continuity);
    end
end
                       
% BUILD THE OUTPUT idmodes
idmodes.s=s;
idmodes.adjacences=adjacences;
names=fieldnames(idmodes1);
% add to idmodes the fields of idmodes1
for i=1:length(names)
    idmodes = setfield(idmodes,names{i},idmodes1.(names{i}));
end


fprintf('\n ***************************************************');
fprintf('\n *** PWA regression terminated successfully !');
fprintf('\n *************************************************** \n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  function hit_compute_modes.m
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [F,idmodes,adjacences]=...
    hit_compute_modes(Xid,yid,ndata,s,xi,inliers,Regressor_set,continuity)

% FIND THE MODE DATASETS F_i 
%
% This is done by using the bijective map between LPVs/FVs
% and datapoints
% -------------------------------------------------------------

global idpar 

F=hit_compute_F(s,inliers,Xid,yid);

% idmodes: structure containing all information on the identified models
%
%	  idmodes.par: cell array: idmodes.par{i} PV of the
%	  i-th mode
%     idmodes.regions: cell array:  idmodes.regions{i} is a cell array.
%                       the region of the i-th mode is
%                       idmodes.regions{i,1}*x \leq idmodes.regions{i,2}
%     idmodes.cov: cell array: idmodes.cov{i} is the covariance of the
%     i-th PV
%     idmodes.regions_sim: cell array defined only if idpar.Regressor_set_sim is
%                           defined. The region for simulating the i-th
%                           mode is
%                       idmodes.regions_sim{i,1}*x \leq idmodes.regions_sim{i,2}

idmodes.par=cell(s,1);
idmodes.regions=cell(s,2);
idmodes.cov=cell(s,1);

% DETERMINE THE SEPARATING HYPERPLANES BETWEEN MODE REGIONS
%
% -------------------------------------------------------------


% sep_hyp: cell array. the ij-entry is the hyperplane [w b]*[x,1]'=0
% separating the cluster i from the cluster j
%

sep_hyp=cell(s,s);

% adjacences: each row is a pair (i,j) indicating that the regions i and j are adjacent
adjacences=[];
% adjacences_hyp: each row contains the vector [w b] of the hyperplane that separates
% region i from region j (the regions are stored in the corresponding row of adjacances)
% according to the formula [w b]*[x,1]'=0
%

adjacences_hyp=[];

if s>1 % process the regions only if there is more than a single mode !
    switch lower(idpar.patt_rec_algo)
        case{'svc'}
            [sep_hyp,idmodes.pattern_rec_valid]=hit_prec_SVC(Xid,F,s);
        case{'mrlp'}
            [sep_hyp,idmodes.pattern_rec_valid]=hit_prec_MRLP(Xid,F,s);
        case{'psvc'}
            [sep_hyp,idmodes.pattern_rec_valid]=hit_prec_PSVC(Xid,F,s);
    end

    %
    %  automatic constraints reduction
    %

    [idmodes.regions,adjacences,adjacences_hyp]=hit_reduce_constraints(sep_hyp,Regressor_set,idpar.LPsolver_cnstr_reduction);
    for i=1:length(idmodes.regions)
        if isempty(extreme(idmodes.regions(i)))
            fprintf('\n ======================================================== ');
            fprintf('\n WARNING: A MODE REGION IS EMPTY');
            fprintf('\n when using idpar.Regressor_set as regressor set.');
            fprintf('\n This has been likely caused by:');
            fprintf('\n 1) A large number of misclassified datapoints');
            fprintf('\n 2) The "true" regions are not polytopic\n');
            fprintf('\n In any case, it is recommended to try a different');
            fprintf('\n pattern recognition algorithm for finding the regions.');
            fprintf('\n ======================================================== \n');
        end
    end
    % idpar.Regressor_set_sim defines the set X on which the identified model will be used
    % this can be bigger than Regressor_set and, if defined, the computation of
    % idmodes.regions_sim is performed
    if isfield(idpar,'Regressor_set_sim') & ~isempty(idpar.Regressor_set_sim)
        fprintf('\n The regressor set for simulation is used for computing the regions.')
        [idmodes.regions_sim,trash1,trash2]=hit_reduce_constraints(sep_hyp,idpar.Regressor_set_sim,idpar.LPsolver_cnstr_reduction);
        for i=1:length(idmodes.regions)
            if isempty(extreme(idmodes.regions(i)))
                fprintf('\n ======================================================== ');
                fprintf('\n WARNING: A MODE REGION IS EMPTY');
                fprintf('\n when using idpar.Regressor_set as regressor set.');
                fprintf('\n This has been likely caused by:');
                fprintf('\n 1) A large number of misclassified datapoints');
                fprintf('\n 2) The "true" regions are not polytopic\n');
                fprintf('\n In any case, it is recommended to try a different');
                fprintf('\n pattern recognition algorithm for finding the regions.');
                fprintf('\n ======================================================== \n');
            end
        end
    end
else % there is only a single mode: set the regions equal to the regressor set
    idmodes.regions=Regressor_set;
    if isfield(idpar,'Regressor_set_sim') & ~isempty(idpar.Regressor_set_sim)
        idmodes.regions_sim=idpar.Regressor_set_sim;
    end
    adjacences=[];
end
%
%  COMPUTE THE MODE PVs
%
% -------------------------------------------------------------

%
%  Find the continuos PWA model if idpar.continuity='c'. If not find the
%  discontinuous PWA model
%
% -------------------------------------------------------------
switch upper(continuity)
    case{'C'}
        [final_PVs_cont]=hit_compute_cont_modes(F,ndata,s,[xi.weights{:}],adjacences,adjacences_hyp);
    case{'D'}
        [final_PVs,idmodes.cov,...
            Weight_primal]=hit_compute_disc_modes(F,ndata,[xi.weights{:}]);
end %ends switch

% stores the coefficients of the models in the structure idmodes

switch upper(continuity)
    case{'C'}
        for i=1:s
            idmodes.par{i}=(final_PVs_cont(i,:))';
        end
    case{'D'}
        for i=1:s
            idmodes.par{i}=(final_PVs(i,:))';
        end
end

%================================================================================
function F=hit_compute_F(s,inliers,Xid,yid)
% F: structure containing mode datasets (i.e. classified datapoints)
% made of inliers
%
%	  F.X: cell array: F.X{i}  matrix of X-points belonging to the i-th
%	  mode
%     F.y: cell array: F.y{i} vector of y-points corresponding to the X-points
%						 belonging to the i-th mode.
%     F.pos: cell array: F.pos{i} vector of indices of the X points in the i-th mode.
%						 Such indices identify the position (i.e. the # of the row)
%						 in the matrix Xid of each point

F.X=cell(s,1);
F.y=cell(s,1);
F.pos=cell(s,1);

% Only the inliers go in the sets F_i
%
% number of inliers
ninl=size(inliers.pos,1);
% cycle over all the inliers and create the mode datasets F
% FROM NOW ON, THE MODE DATASETS CONTAIN ONLY THE INLIERS
for i=1:ninl
    F.X{inliers.class(i)}=[F.X{inliers.class(i)}; Xid(inliers.pos(i),:)];
    F.y{inliers.class(i)}=[F.y{inliers.class(i)}; yid(inliers.pos(i))];
    F.pos{inliers.class(i)}=[F.pos{inliers.class(i)},inliers.pos(i)];
end



