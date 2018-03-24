% filename: initialize_matrixA.m
% Purpose:  creating the boolean matrix of nodes 
% regarding Betweenness centrality or random failures at specific position based on the parametrization
% setting for range 
% Input: 
% - position - x and y coordinates for each network agent
% - param - parametrization struct
% - options - graph options (set the parametrization for graph properties
% computation
% Output: 
% - matrix - weighted adjacency matrix 


function [matrix] = initialize_matrixA(position,param,options)

distance = squareform(pdist(position,'euclidean')); % compute the euclidean distance between each network node

matrix=zeros(size(distance,1),size(distance,1)); 

matrixW = exp(-(distance.^2)/(2*param.sigma^2)); % compute the network weight        

[r,c]= find((param.range) >= triu(distance,0) & triu(distance,0) ~= 0);  % verify which nodes can communicate with the specified range

for i=1:length(r)
    weight = 1;
    if ~options.unweighted
        weight = matrixW(c(i,1),r(i,1)); % path cost
    end

    matrix(r(i,1),c(i,1))=weight;
    matrix(c(i,1),r(i,1))=weight;
end
