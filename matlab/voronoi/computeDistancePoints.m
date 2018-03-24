% filename: computeDistancePoints.m
% Purpose:  verifies spots that are covered by robots considering the distance of each spot to each node into the network, if the distance
% between the spot and the node is lower than the node range, the spot is set as covered by this node
% Input: 
% - XY - spot positions - x and y coordinates
% - position - x and y coordinates for each network agent
% - param - parametrization struct
% Output: 
% matrix A(n,m) where n= number of node and m= number of spots- boolean matrix mapping each node to each spot, with 1 meaning that the spot is covered by the node

function A = computeDistancePoints(XY,position,param)

distance=squareform(pdist([position;XY],'Euclidean'));
D=distance(1:size(position,1),size(position,1)+1:size(distance,2));
covered=find(D<=param.range); %verifies range
A=zeros(size(D,1),size(D,2));
A(covered)=1; 