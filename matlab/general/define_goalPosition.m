% filename: define_goalPosition.m
% Purpose:  for each network, verifies its vulnerability level and the goal position, that is the position the node needs 
% to move toward to mitigate its state of vulnerability.  If a node is not
% vulnerable its goal position is its current position
% Input: 
% - position - x and y coordinates for each network agent
% - param - parametrization struct
% Output: 
% - goal_Position - x and y coordinates for each network agent 
% - vN - struct containing information about the vulnerable nodes  
% Reference :
% C. Ghedini, C. Secchi, C. H. C. Ribeiro, and L. Sabattini. Im-
% proving robustness in multi-robot networks. In Proceedings of
% the IFAC Symposium on Robot Control (SYROCO), Salvador,
% Brazil, aug. 2015.

function [goal_Position vN] = define_goalPosition(position,param)

vN=local_vulnerability(position,param); 

goal_Position=position;

for i=1:size(vN,1)
    
    current_position =vN{i,2};  % current position
    
    goal_Position(vN{i,1},:) = mean(vN{i,3}(:,2:3),1); % new position - average of target nodes   

end
