%generates the benchmark databaset according to the parameter setting -
%setModelDatabase file.

close all;
clear all;

clc;

% model parameter setup

graph_init;   % graph library
graphOptions; % graph parameters
setModelDatabase;    % model parameters

dealCounter=0;

data=[];

while param.network < param.numberNetworks
    
    deal=false;
    
    properties=[];

    position=[];
    
    while dealCounter < param.deals && ~deal
        
        RandStream.setGlobalStream(RandStream('mt19937ar','seed',sum(100*clock)));
         
        [position connected] = graph_generator(param,options);
        
        properties = [compute_properties(position,options,param)];
        
        if (connected && param.connected) && (properties(1,5) >= param.epsilon && param.algebraic)
            deal=true;
        else
            dealCounter = dealCounter + 1;
        end
        
    end 
    if param.algebraic && properties(1,5) < param.epsilon
        disp('It was not possible to generate the defined algebraic connectivity')
        break;
    end
    
    if ~(connected) && param.connected
        disp('It was not possible to generate a connected graph')
        break;
    end
    
    param.network=param.network+1;
    
    disp(sprintf('Network number= %d', param.network));
    
    updated_path=[param.path num2str(param.network) '\' ];
    
    if ~isequal(exist(updated_path, 'dir'),7)
        mkdir(updated_path);
    end
    
    save([updated_path 'position'], 'position');

    save([updated_path 'properties'], 'properties');
    
    data=[data;properties];
        
end
save([param.path 'data'], 'param');
plotPropertiesDatabase(data,param);

free_all;