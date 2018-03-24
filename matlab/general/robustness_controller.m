% filename: robustness_controller.m
% Purpose: computes the weights for improving the robustness to failures
% considering the difference between the node current position and the
% goal position and a linear velocity (param.k)
% - position - x and y coordinates for each network agent
% - goal_Position - goal position for nodes
% - param - parametrization struct
% data: 
%  data - weights for robustness improvement
% Reference:
% C. Ghedini, C. Secchi, C. H. C. Ribeiro, and L. Sabattini. Improving robustness in multi-robot networks. In Proceedings of
% the IFAC Symposium on Robot Control (SYROCO), Salvador, Brazil, aug. 2015.

function [data] = robustness_controller(position,goal_Position,param)

for i=1:size(position,1)
    
     data(i,1) = param.k * (goal_Position(i,1) - position(i,1));
     
     data(i,2) = param.k * (goal_Position(i,2) - position(i,2));
     
end