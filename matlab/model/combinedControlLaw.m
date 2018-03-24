% filename: combinedControlLaw.m
% Purpose: executes the differencial equation for the combined control law
% during the time set
% - position - x and y coordinates for each network agent
% - goal_position - x and y goal position coordinates for each network agent regarding robustness improvement 
% - centroid - x and y goal position coordinates for each network agent
% regarding coverage
% - param - parametrization struct
% - options - graph options (set the parametrization for graph properties
% computation
% Output: 
% factorU - final x and y coordinates for each network agent after
% simulation time interval

function [factorU] = combinedControlLaw(t,positionI,goal_position,centroid,param,options)

persistent data % persist simulation data for supporting the midia generation

nv=round(size(positionI,1)/2);

position=[positionI(1:nv) positionI(nv+1:size(positionI,1))];

[l c] = size(position);

if param.save % persist data
    data=[data;t compute_properties(position,options,param) graph_connected(position,param)];
    save([param.pathR num2str(param.network) '\properties_' num2str(param.idx)], 'data');
end

[A] = initialize_matrixA(position,param,options); % matrix adjacency initialization

ac  =  algebraic_connectivity_New(A,param.normalized); % algebraic connectivity computation

connected=1; %set network as connected

if ~(graph_connected(position,param))  % evaluates if graph is connected  
    if param.save %save iteration data (disconnection) if set
        connectivityData = dir(fullfile([param.pathR num2str(param.network) '\timeD*.*']));        
        if isempty(connectivityData)
            timeDisc=param.idx*t;
            save([param.pathR num2str(param.network) '\timeDisc' num2str(param.idx) '_' num2str(timeDisc)], 'timeDisc');
        end
    end
    param.adpGain=0; %the connectivity gain for this iteration is zero
    connected=0; % flag of not connected graph
end

tolerance = 10E-3; % tolerance for comparing float values 

% CONNECTIVITY MAINTENANCE

% if the algebraic connectivity is lower than epsion its gain is set to
% zero in order to avoid convergence problems. The algebraic connectivity mechanism ensures convergence for ac>epsilon. 
% as networks face perturbation the ac can abruptally decrease lower than
% epsilon. Otherwise, param.adpGain=param.gainConnectivityController

if abs(ac-param.epsilon) < tolerance  
   param.adpGain=0;
else
   param.adpGain=param.gainConnectivityController;
end

if param.adpGain==0 % for evaluation performance improvement
   posxy=zeros(l,c);
else
   [posxy] = connectivity_controller(position,param,options);   % connectivity factor
end

posxy=posxy.*param.adpGain; % connectivity control weight 


% COLLISION AVOIDANCE 

if param.collisonAvoidance && param.gainCollisionAvoidance>0  % if collision avoidance is active and its gain > 0

    dataxy = collisionAvoidance_controller(position,param); % collision avoidance

    dataxy = dataxy.*param.gainCollisionAvoidance;
    
    posxy=[posxy+dataxy];
end

% ROBUSTNESS IMPROVEMENT CONTROL LAW 

dotxy=zeros(l,2); 

if param.gainRobustnessControl~=0 % for evaluation performance improvement 
     
   [dotxy] = robustness_controller(position,goal_position,param);  % robustness factor
 
end

dotxy=dotxy.*param.gainRobustnessControl; % robustness control weight 

% for shifting agents on the cartesian plane

dotxy(:,1) = dotxy(:,1) + (param.constantFormationX * param.gainFormationX);

dotxy(:,2) = dotxy(:,2) + (param.constantFormationY * param.gainFormationY);


%  COVERAGE AREA IMPROVEMENT CONTROL LAW
 

dotxyC=zeros(l,2);

if  param.gainCoverageController>0  % for evaluation performance improvement 
     
   [dotxyC] = coverage_controller(position,centroid,param);   % coverage factor
    
   dotxyC=dotxyC.*param.gainCoverageController; % coverage weight
   
end

xxx = [posxy(:,1) + dotxy(:,1) + dotxyC(:,1); posxy(:,2) + dotxy(:,2) + dotxyC(:,2)]; % final weight computation
 
factorU=xxx;
