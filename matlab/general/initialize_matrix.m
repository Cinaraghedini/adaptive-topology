% filename: initialize_matrix.m
% Purpose:  initialize the network matrices: adjacency, weighted and full weighted (considering fully connected network) for normalization purposes   
% Input: 
% - position - distance matrix 
% - param - parametrization struct
% - options - graph options (param the parametrization for graph properties
% computation
% Output: 
% - matrix - adjacency matrix 
% - matrixW - weighted matrix
% - matrixWN - weighted matrix for a fully connected network


function [matrix matrixW matrixWN] = initialize_matrix(distance,param,options)

matrix=false(size(distance,1));

matrixW=zeros(size(distance,1));

matrixWN = exp(-(distance.^2)/(2*param.sigma^2));   % set the weights according to the algebraic connectiviy definition

matrixWN(logical(eye(size(matrixWN))))=0; % set the weighted matrix considering a fully connected network

[r,c]= find(param.range >= triu(distance,0) & triu(distance,0) ~= 0); % find nodes that are at the minimum distance to create an edge between them

for i=1:length(r)
    matrix(r(i,1),c(i,1))=true;
    matrix(c(i,1),r(i,1))=true;
    weight = 1;
    if ~options.unweighted 
        weight = matrixWN(c(i,1),r(i,1)); 
    end    
    matrixW(r(i,1),c(i,1))=weight;
    matrixW(c(i,1),r(i,1))=weight;
end