% filename: avgDistCentroid.m
% Purpose:  computes the average distance of each node to the network
% centroid and the network radius
% Input: 
% - position - x and y coordinates for each network agent
% - param - parametrization struct
% - options - graph options (set the parametrization for graph properties
% computation
% Output: 
% - R - network radius
% - avgDist - average distance of each node to the network centroid


function [R,avgDist]=avgDistCentroid(position)

[xc1,yc1,R,a] = circfit(position(:,1),position(:,2)); % fits a circle  in x,y plane

distance = squareform(pdist([[xc1,yc1];position],'euclidean'));

avgDist=sum(distance(1,:))/size(position,1);