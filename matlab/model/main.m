% filename: main.m
% Purpose: protocol for evaluating the performance of the combined control strategy:
% algebraic connectivity, robustness to failure and area
% coverage improvement based on a parametrization file (setModel.m)
% Parametrization files :
% - graphOptions.m:  graph parameters
% - setModel.m: model parametrization 
% Output:
% - Files tracking each iteration results for network properties
% - Figures with analitical results

close all;
clear all;
clc;
warning off;

% model parameter setup

graph_init;   % graph library
graphOptions; % graph parameters
init();

setModel;    % model parameters

%network id initialization 

if isempty(param.networkList) % if there is no specific networks set, sequencial numbers considering the number of network are generated
   param.networkList=transpose([1:param.numberNetworks].');
end

param.numberNetworks=length(param.networkList);  % set the number of networks 

auxObstacles=obstacle_generator(param); %generates obstacle positions

param.obstacles=auxObstacles; % set obstacles
        
[param.failureT] = failureTime_generator(param);  %generates failure time distribution

propertiesF =[]; % persists the average of the network evolution for each gain combination 

discT=[]; %persists the average of the disconnection time iteration for each gain combination 

param.legend=[]; % generates the picture legends  

for k=1:size(param.gain,1) % for each gain combination  - it can be used for running specific parameter setting, thus it is possible to run a set of scenarios

    % initialization of iteration variables
    
    discIt=[];
    
    if param.attack % if attacks are performed the number of properties measured match the number of attacks to be perfomed      
       properties=zeros(size(param.failureT,2)+1 ,10);
       coverageData=zeros(size(param.failureT,2)+1 ,18);       
    else  
        properties=zeros(param.iteration+1,10);
        coverageData=zeros(param.iteration+1,18);
    end
    
    % gain setting 
    param.gainConnectivityController=param.gain(k,1);
    param.gainRobustnessControl=param.gain(k,2); 
    param.gainCoverageController=param.gain(k,3);
    
    % set labels for picture and folders/file name generation
    labelP=strcat(num2str(param.gainConnectivityController),'_',num2str(param.gainRobustnessControl),'_',num2str(param.gainCoverageController));
    param.fileName = labelP;     
    param.legend=[param.legend,param.labels(1,k)];
    param.XTickLabel=[param.XTickLabel,{labelP}];
    param.pathR = [param.mainPath labelP '\'];
    
    for i=param.networkList % for each network 
      
        param.listNodeFailure=[]; % persists the nodes that fails
        
        param.network=i; % set this iteration network  id
        
        load([param.path num2str(param.network) '\position']); % load the network i from the database path - position is a matrix containing x and y coordenates for each network agent
     
        % create the folder for results persistence
        
        if ~isequal(exist([param.pathR,num2str(param.network),'\']),7)
            mkdir([param.pathR,num2str(param.network),'\']);
        end
       
        % for tracking the execution status
        disp(sprintf('Network: %d, Connectivity Gain: %d, Robustness Gain: %d, Coverage Gain: %d', i, param.gain(k,1),param.gain(k,2),param.gain(k,3)));
     
        % for perfoming the adaptation process
              
        [propertiesI, dscData, param, coverageDataI] = adapt_network(position,param,options); % adpative mechanism
        
        param.obstacles=auxObstacles; % set obstacles
        
        properties=properties+propertiesI; % for computing the average of the networks properties
        
        coverageData=coverageData+coverageDataI; % for computing the average of the coverage area properties
        
        if ~isempty(dscData) % if the network got disconnected 
            discIt=[discIt;dscData]; %traking disconnection simulation time
        end
        
        % save the i network simulation data for generating videos
        save([param.pathR num2str(param.network) '\properties' ], 'properties');
        save([param.pathR num2str(param.network) '\coverageData'], 'coverageData');
        save([param.pathR num2str(param.network) '\disconnectionI'], 'discIt');
        save([param.pathR num2str(param.network) '\param'], 'param');
        
    end
       
    % this "section" compute and persists average of all network for the specific gain setting  
    
    if ~isequal(exist(param.pathR, 'dir'),7) % create root directory
        mkdir(param.pathR)
    end
    
    if param.numberNetworks>0 
        
        properties=properties/param.numberNetworks; %compute the properties average
        
        coverageData=coverageData/param.numberNetworks; % compute the coverage area average
        
        if ~isempty(discIt) % if some network got disconnected
            discIt=mean(discIt); %take the average disconnection time
        else
            discIt=0;
        end
            
       plotPropertiesEvolution(properties,discIt,param); % plot the average network properties evolution for the specific gain setting
        
       plotCoverageData(coverageData,param); % plot the average coverage area data evolution for the specific gain setting
       
       close all;
        
    end
   
   % save and persist each gain setting simulation data
    
    save([param.pathR 'properties' ], 'properties');
    save([param.pathR 'coverageData' ], 'coverageData');
    save([param.pathR 'disconnectionI'], 'discIt');
    save([param.pathR 'param'], 'param');
        
    propertiesF =[propertiesF properties]; %network properties
    coverageDataF(k,1) ={coverageData}; %coverage data
    discT=[discT discIt]; %disconnection time
    
end

% this section takes care of the total simulation data 

plotProperties_Integrated(propertiesF,discT,param); % plot integrated simulation data for network properties - all the gain settings 
plotCoverage_Integrated(coverageDataF,param);  % plot integrated simulation data for coverage - all the gain settings 

% save integrated data

save([ param.mainPath 'propertiesF'], 'propertiesF');
save([ param.mainPath 'coverageDataF'], 'coverageDataF');
save([ param.mainPath 'discT'], 'discT');
save([ param.mainPath  'param'], 'param');

free_all;