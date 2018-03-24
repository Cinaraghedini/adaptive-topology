% filename: findPos_State.m
% Purpose:  generating a figure from one timestamp data (network state) 
% Input:
% - data: matrix containing the simulation time, node positions, state of node (vulnerable or not) 
% - state: 1 for returning vulnerable nodes, 0 for non-vulnerable
% Output: 
% - dataPlot: matrix with vulnerable/non-vulnerable nodes data 

function [dataPlot] = findPos_State(data,state)

pos=find(data(:,5)==state); 

dataPlot=[];

for i=pos
  dataPlot=[dataPlot; data(i,:)];
end    
