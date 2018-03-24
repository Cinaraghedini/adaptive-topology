% filename: adapt_newtork.m
% Purpose:  for each network and gain setting executes the ODE and compute the necessary information for network adaptation, such as vulnerability level and voronoi cell centroids
% Input:
% - position - x and y coordinates for each network agent
% - param - parametrization struct
% - options - graph options (set the parametrization for graph properties
% computation
% Output:
% - properties - network propertes evolution during the simulation time
% - dscData - disconnection time
% - param - parametrization struct
% - coverageData - area coverage evolution during the simulation time

function [properties, dscData, param ,coverageData] = adapt_network(position,param,options)


%% local variables initialization

properties=[]; %network properties
dscData=[]; % disconnection time
coverageData=[]; % coverage estimation data
sumCoverage=[];
sumProperties=[];

RandStream.setGlobalStream(RandStream('mt19937ar','seed',sum(100*clock))); %random seed generation

idx=param.t0; % initial time setting
labels = num2str([1:size(position,1)]'); % define the agents labels - it is necessary to control de node removed from the network

N=size(position,1); %network size


while idx < param.tf && ~isempty(position) % while there is nodes into the network and the simulation time is not the final one
    
    param.idx=idx; % set the time - is use into the ODE for persisting the iteration results
    
    param.fractionIteration = 1/size(position,1); % fraction of node still into the network - necessary to compute the robustness
    
    propertiesI=compute_properties(position,options,param); %compute the network properties
    
    coveragePropI=coverageProperties(position,param,options); %compute the properties related to the coverage estimation
    
    if param.t0==idx %if is the first iteration - for  keeping the initial value out of the average computation
        properties=[properties ;idx propertiesI]; %compute the network properties
        coverageData=[coverageData; idx coveragePropI];
    end
    
    %%% inserting the network properties and coverage estimation for this iteration into a matrix for analytical results
    if param.attack %the fault-prone scenario is set - (attacks are performed)
        f=1-(size(position,1))/N;   %franction of network in the network
        sumCoverage=[sumCoverage;coveragePropI];
        sumProperties=[sumProperties;propertiesI];
    else %the fault-free scenario is set - (attacks are not performed)
        properties=[properties ;idx propertiesI];
        coverageData=[coverageData; idx coveragePropI];
    end
    
    % information necessary to the robustness controller
    
    if param.gainRobustnessControl>0 %if the robustness control law is active compute for each node its vulnerability and goal position
        [goal_Position vN] = define_goalPosition(position,param); % vN is the list of vulnerable nodes
    else
        goal_Position=position; %if robustness is not active the goal position is set as the current position
        vN=[]; %there is no vulnerable nodes
    end
    
    
    % vArray - list of nodes with 1 for vulnerable nodes
    vArray=zeros(size(position,1),1);
    if ~isempty(vN)
        vArray(cell2mat(vN(:,1)),1)=1;
    end
    
    % information necessary to the coverage controller
    
    if param.gainCoverageController>0
        [centroids]=define_centroid_local(position,options,param); % defines the controid of each node voronoi diagram cell
    else
        centroids=position; % if the coverage gain is not active the current node is position is assumed as the reference for centroid 
    end
    
    % combined control law
    
    clear  combinedControlLaw;
    
    [t, positionI]=ode45(@(t,x) combinedControlLaw(t,x,goal_Position,centroids,param,options),[param.t0 param.ti],position);
    
    position=[transpose(positionI(size(positionI,1),1:round(size(positionI,2)/2))) transpose(positionI(size(positionI,1),round(size(positionI,2)/2)+1:size(positionI,2)))];
    
    if param.save % persists iteration data for videos generation - vulnerable nodes, network properties, node that failed, node labels, node positions
        iterationData=[[-1;[propertiesI(size(propertiesI,1),8)];t] [[transpose(str2num(labels)) transpose(str2num(labels))]; [transpose(vArray) transpose(vArray)]; positionI]];
        save([param.pathR num2str(param.network) '\' param.fileId '_' num2str(idx)], 'iterationData');
    end
    
    
    % for failure-prone scenario
    
    if param.attack && ~isempty(find(param.failureT==idx))
        
        [positionAux pos]=attack_network(position,param,options); % promote perturbantion into the network by removing nodes according the strategy set
        
        if ~isempty(pos)
            
            param.listNodeFailure=[param.listNodeFailure; position(pos,:)]; % update the nodes that failed
            
            param.obstacles=[param.obstacles;position(pos,:)]; %update obstacles - the node which failed is now an obstacle into the environment
            
            labels(pos,:)=[];  % update the node labels
            
        end
        
        properties=[properties ;f mean(sumProperties,1)]; % update the network properties
        coverageData=[coverageData; f mean(sumCoverage,1)];
        sumCoverage=[]; % update the coverage estimation data
        sumProperties=[];
        position=positionAux; % update the node positions, since node(s) failed
        
    end
    
    if ~(graph_connected(position,param)) && isempty(dscData) % update the failure time
        
        dscData=idx*param.ti;
        
    end
    
    idx=idx+param.ti;
    
end
