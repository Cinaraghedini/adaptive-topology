%generates the benchmark databaset according to the parameter setting -
%setModelDatabase file.

close all;
clear all;

clc;

% model parameter setup

graph_init;   % graph library
graphOptions; % graph parameters
setModelObstacles;    % model parameters


while param.network < param.numberNetworks
    
    param.network=param.network+1;
    
    updated_path=[param.path num2str(param.network) '\' ];
    
    if ~isequal(exist(updated_path),7)
        disp(sprintf('Network= %d does not exist', param.network));
    else
        
        disp(sprintf('Network number= %d', param.network));
        
        obstacles = obstacle_generator(param);
                
        save([updated_path 'obstacles'], 'obstacles');
        
    end
end

free_all;