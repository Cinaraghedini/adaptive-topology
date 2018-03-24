% filename: coverageEstimation.m
% Purpose:  computes the spots that are covered, uncovered and served by more than one node (overlapped) by the network nodes 
% Input: 
% A matrix (n,m) where n= number of node and m= number of spots- boolean matrix mapping each node to each spot, with 1 meaning that the spot is covered by the node
% Output: 
% coverageIdx array containing 
% fractionInside  
% - fractionInside(1,1) fraction covered spots 
% - fractionInside(1,2) normalized fraction
% fractionOutside  
% - fractionOutside(1,1) fraction of uncovered spots 
% - fractionOutside(1,2) normalized fraction
% fractionOverlap - fraction of the sum of nodes serving each spot by the sum considering that all spots are served by all nodes.  
 
function [fractionInside, fractionOutside,FractionOverlap] = coverageEstimation(A)

normFactor=1-(size(A,1)/(size(A,1)+size(A,2)));

columnSum=sum(A(:,:));
coveredPoints=length(find(columnSum));
notcoveredPoints=length(find(~columnSum));

fractionInside(1,1)=coveredPoints/size(A,2);
fractionInside(1,2)=fractionInside(1,1)*normFactor;
fractionOutside(1,1)=notcoveredPoints/size(A,2);
fractionOutside(1,2)=fractionOutside(1,1)*normFactor;

FractionOverlap(1,1)=(sum(columnSum)-coveredPoints)/((size(A,2)-1) * size(A,1));