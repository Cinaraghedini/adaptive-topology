function [centers,cost,inl,class,etime,spec] = hit_kmeans(v,opt)
%HIT_KMEANS Weighted KMEANS algorithm.
%
% -------------------------------------------------------------------------
% DESCRIPTION
% -------------------------------------------------------------------------
% [centers,cost,inl,class,etime,spec] = hit_kmeans(v,opt)
%
% -------------------------------------------------------------------------
% INPUT
% -------------------------------------------------------------------------
% v(i,:): i-th point to be clustered. Points are ROW vectors
% opt: structure of specific parameters
%   MANDATORY FIELDS
%   opt.init_centers: 'COVARIANCES'; 'SCALARS'. Specify how to weight
%   points in initializing the centers. 'SCALARS' is fast but 'COVARIANCES'
%   is much more precise.
%   opt.centers: 'COVARIANCES'; 'SCALARS'. Specify how to weight
%   points in updating the centers. 'SCALARS' is fast but 'COVARIANCES'
%   is much more precise.
%   opt.repetitions: n. of times Kmeans must be run. Only the best
%   results will be given as the outputs.
%   opt.s: number of clusters. 
%   DEPENDENT FIELDS
%   opt.w(i): weight for the i-th point.
%   opt.IR{i}: inverse of the covariance associated to the i-th point.
%   OPTIONAL FIELDS
%   opt.options: vector of options for k-means
%   options(1) is set to 1 to display error values; If options(1) is set to
%   0, then only warning messages and cost at each run are displayed.
%   options(1) is -1, then nothing is displayed. Default=-1
%
%   options(2) is a measure of the absolute precision required for the
%   value of centres at the solution.  If the absolute difference between
%   the values of centres between two successive steps is less than
%   options(2), then this condition is satisfied. Default= 1e-4
%
%   options(3) is a measure of the precision required of the error
%   function at the solution.  If the absolute difference between the
%   error functions between two successive steps is less than options(3),
%   then this condition is satisfied. Both this and the previous
%   condition must be satisfied for termination. Default= 1e-4
%
%   options(4) is the maximum number of iterations in a single run of
%   Kmeans; Default= 100.
%
% -------------------------------------------------------------------------
% OUTPUT
% -------------------------------------------------------------------------
% centers(i,:): center of the i-th cluster (is a ROW vector) found in
% the best run.
%
% cost: value of the clustering cost functional in the best run.
%
% inl: indexes of datapoints that are inliers after clustering (i.e.
% not discarded by the clustering algorithm). Since Kmeans oes not perform
% outlier detection, all input points will be inliers.
%
% class(i): classification of the i-th inlier. class(i)=j means
% that the i-th inlier belongs to the j-th cluster.
%
% etime: elapsed time for the execution of ALL the Kmeans runs.
%
% spec: structure with special outputs
%   spec.costs is the vector of length opt.repetitions storing the
%   costs at each iteration of Kmeans in the best run of the algorithm.
%   spec.cost_each_clust: vector of costs for each cluster found (in the
%   best run).
% -------------------------------------------------------------------------
% ACKNOWLEDGMENTS
% -------------------------------------------------------------------------
% This unction is based on the kmeans.m routines in the NetLab toolbox
% developed at the Aston University (Birmingham).


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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% CHECK FIELDS OF opt
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check mandatory fields of opt
M={{'init_centers'},{'centers'},{'s'},{'repetitions'}};
for i=1:length(M)
    if ~isfield(opt,M{i})
        error('hit_kmeans: ERROR: one mandatory field of opt is missing')
    end
end

% Check dependent fields of opt
% each row of D={D1;D2;...;Dn} is made of  a triplets of cells
% where D(i,:)={
% {'mandatory field1','mandatory field2'...},
% 'value mandatory field',
% 'dependent field'
%}
% This means that if 'mandatory field1' or 'mandatory field2'
% has the value 'value mandatory field''
% then there must exist 'dependent field'

D={
    {'init_centers','centers'},'covariances','IR'
    {'init_centers','centers'},'scalars','w'
    };
% # rows of D
rD=size(D,1);
for i=1:rD
    % number of field names in D{i,1}
    lrD=size(D{i,1},2);
    for j=1:lrD
        %field name to compare
        toch=D{i,1}{j};
        %check if the madatory field has the value D{i,2}
        if strcmp(lower(getfield(opt,toch)),lower(D{i,2}))
            % check if the dependent field exists
            if ~isfield(opt,D{i,3})
                error('hit_kmeans: ERROR: a dependent field of opt is missing')
            end
        end
    end
end

% Set optional fields to their default values if not declared
% each row of O={O1;O2;...;On} is made of  a triplets of cells
% where O(i,:)={
% {'optional field i',default value},
%}
% This means that if 'optional field i' does not exist
% it is set to the default value
O={
    {'options',[-1,1e-4,1e-4,100]}
    };
% # rows of O
rO=size(O,1);
for i=1:rO
    if ~isfield(opt,O{i}{1})
        opt=setfield(opt,O{i}{1},O{i}{2});
        % fprintf('hit_kmeans: Field opt.%s set to default value\n',O{i}{1})
    end
end

clear M D O rD lrD rO toch i j

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% END OF CHECKING FIELDS OF opt
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% sizes of the points
[nv, dim] = size(v);

% initialize the vector of the costs
spec.costs=zeros(opt.repetitions,1);

old_e=inf;
if ~all(opt.s <= nv)
    error('hit_kmeans: ERROR: More centers than data')
end

besterror=inf;

%start counter for etime
tic;
fprintf('\n Executing Kmeans with %d clusters...',opt.s);
for lll=1:opt.repetitions

    % initialize the centers
    % sopt1: options for ini_centers

    % number of centers
    sopt1.s=opt.s;
    switch lower(opt.init_centers)
        case{'scalars'}
            sopt1.init_centers='scalars';
            sopt1.w=opt.w;
        case{'covariances'}
            sopt1.init_centers='covariances';
            sopt1.IR=opt.IR;
    end
    centers=ini_centers(v,sopt1);
    if ~all(all(isfinite(centers)))
        fprintf('\n ======================================================== ');
        fprintf('\n ERROR WHILE INITIALIZING THE CENTERS IN KMEANS:');
        fprintf('\n not all entries of the centers are finite.');
        fprintf('\n A possible cause is that scalar weights are "too big"');
        fprintf('\n or confidence matrices are "too small".');
        fprintf('\n (depending on which confidence measure you have chosen).');
        fprintf('\n Both phenomena are likely to happen because of almost');
        fprintf('\n noiseless data.');
        fprintf('\n Rough remedy: add a little noise to output data.');
        fprintf('\n hit_regression.m ends here. ');
        fprintf('\n ======================================================== \n');
        error('hit_regression:end','hit_regression ends here.')
    end
    clear sopt1

    % Sort out the opt.options
    if (opt.options(4))
        niters = opt.options(4);
    else
        niters = 100;
    end
    % set the options for dist2_mod
    switch lower(opt.centers)
        case{'scalars'}
            sopt.centers='scalars';
            sopt.w=opt.w;
        case{'covariances'}
            sopt.centers='covariances';
            sopt.IR=opt.IR;
    end
    % store errors at each iteration (for a fixed sss and in a given
    % repetition)
    errlog = zeros(1, niters);
    % Matrix to make unit vectors easy to construct
    id = eye(opt.s);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % MAIN LOOP FOR A SINGLE RUN OF KMEANS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    exit_flag=0;
    % iteration number
    n=1;
    while(~exit_flag)
        % Save old centres to check for termination
        old_centers = centers;
        % Calculate distances based on existing centers
        d2 = dist2_mod(v, centers,sopt);
        % Assign each point to nearest center
        [minvals, index] = min(d2', [], 1);
        post = id(index,:);
        num_points = sum(post, 1);
        % Adjust the centers based on new distances
        for j = 1:opt.s
            if (num_points(j) > 0)
                temp=find(post(:,j));
                switch lower(opt.centers)
                    case{'scalars'}
                        centers(j,:)=sum(v(temp,:).*(opt.w(temp)*ones(1,dim)), 1)./sum(opt.w(temp));
                    case{'covariances'}
                        AA=zeros(dim,dim);
                        bb=zeros(dim,1);
                        for kk=1:length(temp)
                            AA=AA+opt.IR{temp(kk)};
                            bb=bb+opt.IR{temp(kk)}*(v(temp(kk),:))';
                        end
                        centers(j,:) = (AA\bb)';
                end %ends switch
            end %ends if
        end %ends for j
        % check if centers entries are finite numbers
        if ~all(all(isfinite(centers)))
            fprintf('\n ======================================================== ');
            fprintf('\n ERROR WHILE UPDATING THE CENTERS IN KMEANS:');
            fprintf('\n not all entries of the centers are finite.');
            fprintf('\n A possible cause is that scalar weights are "too big"');
            fprintf('\n or confidence matrices are "too small".');
            fprintf('\n (depending on which confidence measure you have chosen).');
            fprintf('\n Both phenomena are likely to happen because of almost');
            fprintf('\n noiseless data.');
            fprintf('\n Rough remedy: add a little noise to output data.');
            fprintf('\n hit_regression.m ends here. ');
            fprintf('\n ======================================================== \n');
            error('hit_regression:end','hit_regression ends here.')
        end
        % Error value is the total squared distance from cluster centres
        % (sopt still stores the options for calling dist2_mod in
        % the correct way)
        d2 = dist2_mod(v,centers,sopt);
        % Assign each point to nearest centre
        [minvals, index] = min(d2', [], 1);
        e_each_clust=[];
        for kk=1:opt.s
            tm=find(index==kk);
            % error of each cluster
            e_each_clust=[e_each_clust;sum(minvals(tm))];
        end
        % total error
        e = sum(minvals);
        errlog(n) = e;
        if opt.options(1) > 0
            fprintf(1, 'N. centers %4d Cycle %4d  Error %11.6f\n',opt.s(sss), n, e);
        end
        % Test for termination
        if max(max(abs(centers - old_centers))) < opt.options(2) & ...
                abs(old_e - e) < opt.options(3)
            exit_flag=1;
        end
        if n==niters
            if (opt.options(1) >= 0)
                disp('hit_kmeans: Warning: Maximum number of iterations has been exceeded.');
            end
            exit_flag=1;
        end
        % increment the cycle number
        n=n+1;
        old_e = e;
    end % ends while loop
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % END OF A SINGLE RUN OF KMEANS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Store the cost of run n. lll with opt.s centers
    spec.costs(lll,1)=e;
    if e<besterror %take the best vectors
        bestcenters=centers;
        bestpost=post;
        besterror=e;
        beste_each_clust=e_each_clust;
    end
end % ends for lll=1:opt.repetitions


% set the best cost after all the iterations of kmeans for the same
% number of centers
cost=besterror;
if ~isfinite(besterror)
    fprintf('\n ======================================================== ');
    fprintf('\n ERROR WHILE RUNNING KMEANS:');
    fprintf('\n the best error is not finite.');
    fprintf('\n A possible cause is that scalar weights are "too big"');
    fprintf('\n or confidence matrices are "too small".');
    fprintf('\n (depending on which confidence measure you have chosen).');
    fprintf('\n Both phenomena are likely to happen because of almost');
    fprintf('\n noiseless data.');
    fprintf('\n hit_regression.m ends here. ');
    fprintf('\n ======================================================== \n');
    error('hit_regression:end','hit_regression ends here.')
end
% STORE THE BEST CLASSIFICATION WITH opt.s(sss) CENTERS
% store the cost of each cluster
spec.cost_each_clust=beste_each_clust;
% store the classification of vectors for a given number of cluster sss
class_sss{1}=[];
for i=1:nv
    class_sss{1}=[class_sss{1}; find(bestpost(i,:)>0.5)];
end
% store the best centers
centers_sss{1}=bestcenters;
if opt.options(1)==0
    fprintf('Clustering cost at each run for N. centers = %4d:\n',opt.s)
    fprintf('%f\n',spec.costs);
end
% check if there is some empty cluster: if yes, set all silhouette
% coefficients to -1
flag_empty_cluster=0;
for i=1:opt.s
    if isempty(find(class_sss{1}==i));
        flag_empty_cluster=1;
    end
end
if flag_empty_cluster
    fprintf('\n ======================================================== ');
    fprintf('\n WARNING AFTER RUNNING KMEANS: some clusters are empty.');
    fprintf('\n (they will be automatically removed).');
    fprintf('\n The number of clusters that kmeans has to find');
    fprintf('\n is likely to be too high.');
    fprintf('\n ======================================================== \n');
end


% build output arguments
inl=[1:nv]';

% here outputs are formatted
class=[class_sss{1}];
centers=[centers_sss{1}];
etime=toc;
%==========================================================================
function [init_cent]=ini_centers(v,sopt)
%==========================================================================
%
% Title:       ini_centers.m
%
% Project:     Identification of PWA maps
%
% Purpose:     initialize the clusters in the feature space according
%              to the algorithm specified in the string algo
%
%
% Input arguments:
%		v: Matrix of points used for choosing the centers
%		sopt: structure of options
%
% Ouput arguments:
% 		init_cent: Matrix of estimated centers. Each row is a center
%
% Comments:
%
%===============================================================================

% extract the number and size of feature vectors
ndata=size(v,1);
dim = size(v,2);

init_cent=[];
cluster_init={};

index_init=randperm(ndata); % randomize the indices of the feature vectors
v_random=v(index_init,:); % randomize the the feature vectors
size_init=floor(ndata/sopt.s); % assign approximately 1/sopt.s feature vectors to each cluster

switch lower(sopt.init_centers)
    case {'scalars'}
        w_random=sopt.w(index_init);
        % use the weights to compute the centers of the clusters
        for i=1:sopt.s-1
            indices=(i-1)*size_init+1:i*size_init;
            cluster_init{i}=v_random(indices,:);
            init_cent(i,:) = sum((w_random(indices)*...
                ones(1,dim)).*cluster_init{i}, 1)/sum(w_random(indices));
        end

        % the last clusters is created in a separate way in order to collect all
        % the remaining datapoints
        indices=(sopt.s-1)*size_init+1:ndata;
        cluster_init{sopt.s}=v_random(indices,:);
        init_cent(sopt.s,:) = sum((w_random(indices)*...
            ones(1,dim)).*cluster_init{sopt.s}, 1)/sum(sopt.w(indices));

    case {'covariances'}
        % use the Linear System to compute the centers of the clusters

        % randomize in the same way the variances of the
        % feature vectors

        IR_random={sopt.IR{index_init}};

        for i=1:sopt.s-1
            indices=(i-1)*size_init+1:i*size_init;
            cluster_init{i}=v_random(indices,:);
            AA=zeros(dim,dim);
            bb=zeros(dim,1);
            % cycle over all the vectors in cluster_init{i}
            % kk are the absolute indices (the corresponding index
            % in cluster_init{i} is given by kk-(i-1)*size_init )
            for kk=indices
                AA=AA+IR_random{kk};
                bb=bb+IR_random{kk}*cluster_init{i}(kk-(i-1)*size_init,:)';
            end
            init_cent(i,:) = (AA\bb)';
        end

        % the last clusters is created in a separate way in order to collect all
        % the remaining datapoints
        indices=(sopt.s-1)*size_init+1:ndata;
        cluster_init{sopt.s}=v_random(indices,:);
        AA=zeros(dim,dim);
        bb=zeros(dim,1);
        % cycle over all the vectors in cluster_init{sopt.s}
        % kk are the absolute indices (the corresponding index
        % in cluster_init{sopt.s} is given by kk-(sopt.s-1)*size_init )
        for kk=indices
            AA=AA+IR_random{kk};
            bb=bb+IR_random{kk}*cluster_init{sopt.s}(kk-(sopt.s-1)*size_init,:)';
        end
        init_cent(sopt.s,:) = (AA\bb)';
end% ends switch

function n2 = dist2_mod(x, c, sopt)

%===============================================================================
% Title:       dist2_mod.m
%
% Project:     Identification of PWA maps
%
% Purpose:     evaluate the weighted squared norm between two sets of points
%
%
%
% Syntax n2 = dist2_mod(x, c, sopt)
%
% Input arguments:
%		x: first Matrix of points. Each point is a row
%		c: second Matrix of points. Each point is a row
%		opt: structure of specific parameters
%            MANDATORY FIELDS
%            opt.centers: 'COVARIANCES'; 'SCALARS'
%            DEPENDENT FIELDS
%            opt.w: weights
%            opt.IR: inverse of covariances
%            OPTIONAL FIELDS
%            (none)
%
% Output arguments: n2(i,j)=(x_i - c_j) ' * V_{x_i} * (x_i- c_j)
%
% Comments: This function is based on the dist2.m function in the NetLab toolbox
%			developed at the Aston University (Birmingham).
%
%===============================================================================

[ndata, dimx] = size(x);
[ncentres, dimc] = size(c);
if dimx ~= dimc
    error('Data dimension does not match dimension of centres')
end
n2=zeros(ndata,ncentres);
switch(lower(sopt.centers))
    case{'scalars'}
        for i=1:ndata
            V=sopt.w(i);
            dato=x(i,:);
            for j=1:ncentres
                centro=c(j,:);
                temp=dato-centro;
                n2(i,j)=temp*V*temp';
            end
        end
    case{'covariances'}
        for i=1:ndata
            V=sopt.IR{i};
            dato=x(i,:);
            for j=1:ncentres
                centro=c(j,:);
                temp=dato-centro;
                n2(i,j)=temp*V*temp';
            end
        end
end
