% filename: clusterCoefficientCentrality.m
% Purpose:  computes the cluster coefficient of a network node 
% Input: 
% - g - network graph
% - v - node 
% Output: 
% - C -  node v cluster coefficient

function C = clusterCoefficientCentrality(g,v)

% Calculate Cluster Coefficient Centrality
C = 0;
listN = neighbors(g,v); % neighbors of v
[sizeL sizeC]=size(listN);
if (sizeC > 1) % if there is more than neighbor
    neighbor=0;
    for u = 1:sizeC 
        listNn= neighbors(g,listN(u)); % neighborhood of neighbor u of v 
        [linha coluna] = size(listNn);
         for i = 1:coluna
             found=find(listN == listNn(i));
             if(size(found) > 0) % if v and u have i as neighbor
                neighbor = neighbor+1; 
             end
         end
    end
    C  = neighbor/(sizeC*(sizeC-1)); % normalized cluster coefficient
end