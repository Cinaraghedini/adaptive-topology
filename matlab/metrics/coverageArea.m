% filename: coverageArea.m
% Purpose:  computes the area coverage estimation 
% Input: 
% - position - x and y coordinates for each network agent
% - param - parametrization struct
% Output: 
% - area - coverage area rate
% Reference
%Saeed Sedighian Kashi and Mohsen Shari?. Coverage rate calculation in wireless sensor networks. Computing, 94(11):833–856,2012.


function area  = coverageArea(position,param)    

robots = [position repmat(param.range, length(position), 1)];

[area, timeInSeconds] = calcArea(robots);
