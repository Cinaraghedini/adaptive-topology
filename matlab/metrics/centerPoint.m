% filename: centerPoint.m
% Purpose:  given the node positions, creates a square area and uses the
% circfit function for fitting a circle and find the network centroid and radius % 
% centroid and the network radius
% Input: 
% - position - x and y coordinates for each network agent
% computation
% Output: 
% - xc - centroid x position 
% - xy - centroid y position
% - R - radius 
% - a - coeficient describing the circle's equation

function [xc, yc, R]= centerPoint(position)

xlimUpper=max(position(:,1));
xlimLower=min(position(:,1));

ylimUpper=max(position(:,2));
ylimLower=min(position(:,2));

xy=[xlimLower ylimLower; xlimLower ylimUpper; xlimUpper ylimUpper; xlimUpper ylimLower]; 

[xc,yc,R,a] = circfit(xy(:,1),xy(:,2));