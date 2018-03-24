% filename: eLocal.m
% Purpose:  computes the local efficiency of a graph which is the average
% of each node local efficiency
% Input: 
% - mat - adjacency matrix
% - mwI - distance matrix
% - Options - graph options
% Output: 
% - e -  network local efficiency
% Reference :
% Phys Rev Lett. 2001 Nov 5;87(19):198701. Epub 2001 Oct 17.
% Efficient behavior of small-world networks.
% Latora V1, Marchiori M.


function eLocal = eLocal(mat,mw,options)

eLocal = sum(eLocalN(mat,mw,options)) * (1/size(mat,1)); %average of each node local efficiency
   