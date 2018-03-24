% filename: coverage_controller.m
% Purpose: computes the weights for improving the coverage area 
% considering the difference between the node current position and the
% goal position and a linear velocity (param.k)
% - position - x and y coordinates for each network agent
% - goal_Position - goal position for nodes - the centroids of the node voronoi
% cell
% - param - parametrization struct
% data: 
%  data - weights for coverage area improvement


function [dotxy] = coverage_controller(position,goal_Position,param)

for i=1:size(position,1)
    
     dotxy(i,1) = param.k * (goal_Position(i,1) - position(i,1));
     
     dotxy(i,2) = param.k * (goal_Position(i,2) - position(i,2));
     
end
