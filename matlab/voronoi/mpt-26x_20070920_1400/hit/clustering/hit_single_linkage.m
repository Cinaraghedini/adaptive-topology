function [centers,cost,inl,class,etime,spec] = hit_single_linkage(V,OPT)
%HIT_SINGLE_LINKAGE Single-linkage clustering algorithm.
%
% -------------------------------------------------------------------------
% DESCRIPTION
% -------------------------------------------------------------------------
% [centers,cost,inl,class,etime,spec] = hit_single_linkage(V,OPT)
%
% -------------------------------------------------------------------------
% INPUT
% -------------------------------------------------------------------------
% V(i,:): i-th point to be clustered. Points are ROW vectors
% opt: structure of specific parameters
%   MANDATORY FIELDS
%   opt.guess: guessed minimal distance between clusters
%   OPTIONAL FIELDS
%   opt.plot_steps: 'Y' 'N' (set by default to 'N' if absent). It specify
%   if the algorithm has to stop at each iteration and plot the clusters
%   found. The algorithm will continue after that a key is pressed.
%
% -------------------------------------------------------------------------
% OUTPUT
% -------------------------------------------------------------------------
% centers(i,:): center of the i-th cluster (is a ROW vector).
% cost: value of the clustering cost functional.
% inl: indexes of datapoints that are inliers after clustering (i.e.
% not discarded by the clustering algorithm). Note that single-linkage
% performs outlier detection.
% class(i): classification of the i-th inlier. class(i)=j means
% that the i-th inlier belongs to the j-th cluster.
% etime: elapsed time for the execution of single-linkage.
% spec: structure with special outputs
%   spec.costs(i) is the costs of the i-th cluster found.

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

global v d opt C plotpar

v=V;
opt=OPT;
clear V OPT

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% CHECK FIELDS OF opt
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check mandatory fields of opt
M={{'guess'}};
for i=1:length(M)
    if ~isfield(opt,M{i})
        error('hit_single_linkage: ERROR: one mandatory field of opt is missing')
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
    {'plot_steps','N'}};
% # rows of O
rO=size(O,1);
for i=1:rO
    if ~isfield(opt,O{i}{1})
        opt=setfield(opt,O{i}{1},O{i}{2});
        % fprintf('hit_single_linkage: Field opt.%s set to default value\n',O{i}{1})
    end
end

clear M D O rD lrD rO toch i j

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% END OF CHECKING FIELDS OF opt
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%start counter for etime
tic;

% Compute all the distances between points
[nv,dim]=size(v);

%indices of the points
pos=[1:nv]';

% initialize the global variable d
all_dist2(pos)

% Initialize the clusters
% C{i,1}= matrix of points in the i-th cluster. Ex: C{i,1}=[3 7 4]'  means that
%         v(C{i,1},:) are the points in C_1
% C{i,2}= vector of distances between the i-th and the successive clusters
%         C{i,2}(j)= k means that d(k,:) is the distance between the i th and i+j th clusters
%         in fact the distance between two clusters i, i+j is the min of the distances between
%         the points in C{i,1} and the points in C{i+j,1}
%

C=cell(nv,2);

st=1;
for i=1:nv-1
    %each cluster contains one point
    C{i,1}=i;
    en=st+nv-1-i;
    C{i,2}=[st:en]';
    st=en+1;
end
%last cluster
C{nv,1}=nv;
C{nv,2}=[];

%initialize number of clusters
nc=nv;


%create list dc of distances between clusters for the first time
% (after, is updated at the end of the while loop)
%each row is [distance(C_i,C_j),i,j] for j>i

dc=[];
for i=1:nc-1
    distc=d(C{i,2},1);
    indc=[ones(nc-i,1)*i,[i+1:nc]'];
    dc=[dc;distc,indc];
end
% sort dc and take only the first distance
[temp, i_sorted]=sort(dc(:,1));
dc=dc(i_sorted(1),:);
while (~isempty(dc) & dc(1)<=opt.guess)
    %decide the two clusters mi mj to merge
    mi=dc(2);
    mj=dc(3);

    % merge the clusters

    mergec(mi,mj);
    % update nc
    nc=nc-1;
    switch lower(opt.plot_steps)
        case{'y'}
            % START PLOTTING COMMANDS
            %%%%%%%%%%%%%%%%%%%%%%%%%
            tmpcenters=[];
            nclust=nc;
            % Compute centers
            for i=1:nc
                center=[];
                npoints=length(C{i,1});
                AA=zeros(dim,dim);
                bb=zeros(dim,1);

                AA=length(C{i,1});
                bb=v(C{i,1},:);
                center=sum(bb,1)/AA;
                tmpcenters=[tmpcenters;center];
            end

            ninl=size(v,1);
            inl=zeros(ninl,1);

            class=zeros(ninl,1);

            % j: current index of the Mcluster
            jjj=1;
            for i=1:nclust
                % extract the indices of the FVs in Inl and belonging to the i-th Mcluster
                tpos=C{i,1};
                inl(tpos)=tpos;
                class(tpos)=i;
            end
            hit_plot_clust(v,tmpcenters,inl,class,plotpar.marker_modes,opt.plot_fig)
            pause
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % END PLOTTING COMMAND

    % Compute the stop criterion for a cluster

    %create list dc of distances between clusters
    %each row is [distance(C_i,C_j),i,j] for j>i
    %olddc=dc(1);
    dc=[];
    for i=1:nc-1
        distc=d(C{i,2},1);
        indc=[ones(nc-i,1)*i,[i+1:nc]'];
        dc=[dc;distc,indc];
    end
    % sort dc and take only the first distance
    if ~isempty(dc)
        [temp, i_sorted]=sort(dc(:,1));
        dc=dc(i_sorted(1),:);
    end
end

% create output arguments


centers=[];
costs=[];
% Compute centers and cost
for i=1:nc
    center=[];
    npoints=length(C{i,1});
    AA=zeros(dim,dim);
    bb=zeros(dim,1);

    AA=length(C{i,1});
    bb=v(C{i,1},:);
    center=sum(bb,1)/AA;
    centers=[centers;center];

    % Compute the cost of the cluster
    mcost=0;
    for kk=1:npoints
        mcost=mcost+(center-v(C{i,1}(kk),:))*(center-v(C{i,1}(kk),:))';
    end
    costs=[costs;mcost];
end
cost=sum(costs);
spec.costs=costs;
% all vectors are inliers
inl=[1:nv]';
class=zeros(nv,1);
for i=1:nc
    indices=C{i,1};
    class(indices)=i;
end
etime=toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  function mergec.m
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function mergec(mi,mj)

global C d

% let mi<mj so that the merged cluster will have index mi
if mi>mj
    temp=mi;
    mi=mj;
    mj=temp;
end

% number of clusters
nc=size(C,1);

% compute the indices of the v in the merged cluster
mindicesv=[C{mi,1};C{mj,1}];

% extract the distances from C_i
% row kk of pdCi is [index ,kk] meaning that the distance
% between C_i and C_kk is given by d(index,1)
pdCi=[];
% extract the distances between C_kk and C_i, kk=1,...,i-1
for kk=1:mi-1
    pos=C{kk,2}(mi-kk);
    pdCi=[pdCi; pos,kk];
end
% append the distances between C_kk and C_i, kk=i+1..nc
pdCi=[pdCi;C{mi,2},[mi+1:nc]'];

% extract the distances from C_j
% row kk of pdCj is [index ,kk] meaning that the distance
% between C_j and C_kk is given by d(index,1)
pdCj=[];
% extract the distances between C_kk and C_j, kk=1,...,j-1
for kk=1:mj-1
    pos=C{kk,2}(mj-kk);
    pdCj=[pdCj; pos,kk];
end
% append the distances between C_kk and C_j, kk=j+1..nc
pdCj=[pdCj;C{mj,2},[mj+1:nc]'];

% compute the distances between the merged clusters and the others
% each row of pmdC is [index_d,kk] meaning that the distance
% between the merged cluster and C_kk is given by d(index_d,1)

pmdC=[];
for kk=1:nc
    if kk~=mi & kk~=mj
        posdCiCkk=pdCi(find(pdCi(:,2)==kk),1);
        posdCjCkk=pdCj(find(pdCj(:,2)==kk),1);
        [temp,index]=min([d(posdCiCkk,1),d(posdCjCkk,1)]);
        if index==1
            pmdC=[pmdC;posdCiCkk,kk];
        else
            pmdC=[pmdC;posdCjCkk,kk];
        end
    end
end

% STEP1: Update the distances between C_kk and C_i, kk=1,...,i-1;
%        Cut also the distances between C_kk and C_j
for kk=1:mi-1
    ind_upd=find(pmdC(:,2)==kk);
    % upd store the new distance between C_kk and the merged cluster
    upd=pmdC(ind_upd,1);
    C{kk,2}(mi-kk)=upd;
    % cut the distance with the cluster mj
    topreserve=setdiff([1:length(C{kk,2})],mj-kk);
    C{kk,2}=C{kk,2}(topreserve);
end

% STEP2: Cut the distances between C_kk and C_j for kk=i+1,...j-1
for kk=mi+1:mj-1
    % cut the distance with the cluster mj
    topreserve=setdiff([1:length(C{kk,2})],mj-kk);
    C{kk,2}=C{kk,2}(topreserve);
end

% STEP3: Update C_i in order to store the new cluster
C{mi,1}=mindicesv;
if ~isempty(pmdC)
    [indices]=find(pmdC(:,2)>mi);
    C{mi,2}=pmdC(indices,1);
end

% STEP4: Copy C_(kk+1) in C_kk for kk=j,...,nc-1
for kk=mj:nc-1
    C{kk,1}=C{kk+1,1};
    C{kk,2}=C{kk+1,2};
end

% STEP5: Cut C_nc
C1=cell(nc-1,2);
for kk=1:nc-1
    C1{kk,1}=C{kk,1};
    C1{kk,2}=C{kk,2};
end
C=C1;
clear C1


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  function all_dist2.m
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function all_dist2(pos)

global v d opt


[nv,dim]=size(v);

% construct the vector of distances
% dist(k,:)=[||xi_i-xi_j||, pos(i), pos(j)], j>i

d=[];
for j=1:nv
    for k=j+1:nv
        diffjk=v(j,:)-v(k,:);
        d=[d;[norm(diffjk,2),j,k]];
    end
end



