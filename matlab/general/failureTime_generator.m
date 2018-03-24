% filename: failureTime_generator.m
% Purpose:  Generates discrete random failure times 
% Input: parametrization struct - the number of failures
% (param.nrFailures), the initial (param.t0) and the final (param.tf) time
% interval
% Output: failureT - array with failure times.   

function [failureT] = failureTime_generator(param) 

while true
    failureT = sort(randi([param.t0,param.tf],[1,param.nrFailures]));
    if length(unique(failureT)) ==  length(failureT)  %verify if there are repetead numbers
        break;
    end
end
   