% filename: initialize_matrixAdj.m
% Purpose:  creating the boolean matrix of nodes 
% regarding Betweenness centrality or random failures at specific position based on the parametrization
% setting for range 
% Input: 
% - position - x and y coordinates for each network agent
% - param - parametrization struct
% Output: 
% - matrix - boolean adjacency matrix 

function [matrix] = initialize_matrixAdj(position,param)

matrix=false(size(position,1));

distance = squareform(pdist(position,'euclidean')); % compute the euclidean distance between each network node

[r,c]= find(param.range >= triu(distance,0) & triu(distance,0) ~= 0);  % verify which nodes can communicate with the specified range

for i=1:length(r)
    matrix(r(i,1),c(i,1))=true;
    matrix(c(i,1),r(i,1))=true;
end
