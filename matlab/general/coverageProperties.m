% filename: coverageProperties.m
% Purpose:  computes the area coverage estimation 
% Input: 
% - position - x and y coordinates for each network agent
% - param - parametrization struct
% - options - graph options (set the parametrization for graph properties
% computation
% Output: 
% - properties - matrix (1,m) where m is the number of properties evaluated
%  properties(1,1) = network density
%  properties(1,2) = coverage rate
%  properties(1,3) = average coverage rate
%  properties(1,4) = average distance to centroid and Radius
%  properties(1,5:17)= covered, uncovered spots and overlapping areas - see coverageEvaluation for details


function properties=coverageProperties(position,param,options)

n=size(position,1); 

distance = squareform(pdist(position,'euclidean')); % computes the euclidean distante between each pair of of nodes

[m mw mwI] = initialize_matrix(distance,param,options);

properties(1,1)=link_density(m); % network density
properties(1,2)=coverageArea(position,param);  %coverage rate
properties(1,3)=properties(1,2)/n;
[R,avgDist]=avgDistCentroid(position); % average distance to centroid and Radius
properties(1,4)=avgDist; 
properties(1,5:17)= coverageEvaluation(position,param,options); % see coverageEvaluation for details