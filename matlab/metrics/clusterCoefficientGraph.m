% filename: clusterCoefficientGraph.m
% Purpose:  computes the cluster coefficient of a graph which is the average
% of each node cluster coefficient
% Input: 
% - g - network graph 
% Output: 
% - C -  network cluster coefficient

function C = clusterCoefficientGraph(g)

[nv ne] = size(g);
C=0;
if ne>0 %there is more than one edge into the network 
    for v=1:nv % for each network node
        C = clusterCoefficientCentrality(g,v) + C; % node v cluster efficiency
    end
    C = C/nv; % Cluster Coefficient Global - average 
end