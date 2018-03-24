% filename: obstacle_generator.m
% Purpose:  Generates random obstacle positions in a bounded cartesian area
% Input: parametrization struct - the area boundaries (param.lowerX, param.upperX, param.lowerY,
% param.upperY) and the number of obstacle to be generated (param.obstacleN) 
% Output: obstacles - a (n,2) matrix containing x (col 1) and y (col 2) representing the position of each of the n obstacles.  

function [obstacles] = obstacle_generator(param) 

obstacles = [((param.lowerX - param.upperX).*rand(param.obstacleN,1) + param.upperX)  ((param.lowerY - param.upperY).*rand(param.obstacleN,1) + param.upperY)];
