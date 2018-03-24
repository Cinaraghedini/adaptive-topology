% filename: robustness.m
% Purpose:  computes the newtork robustness to failure of central elements
% regarding Betweenness centrality or random failures
% Input: 
% - position - x and y coordinates for each network agent
% - param - parametrization struct
% - options - graph options (set the parametrization for graph properties
% computation
% Output: 
% - data - array of size 1,2 
% - data(1,1) = network robustness
% - data(1,2) = giant component of the remaining network
% Improving robustness in multi-robot networks
% C Ghedini, C Secchi, CHC Ribeiro, L Sabattini
% Proc. 11th IFAC Symp. Robot Control SYROCO, 63-68

function [data]=robustness(position,param,options)

nd=(size(position,1)); %number of nodes into the network

nN=nd;

nrNodes=round(param.fractionIteration*nd); % number of node to be removed at each iteration

data=zeros(1,2);

if nrNodes<nd % if there is enough node to be removed
    
    count=0;
    
    g=graph; % create a graph
    
    connected=true; % consider network as connected
    
    while connected && nd>nrNodes
        
        matrixW = initialize_matrixAdj(position,param); % create the adjacency matrix based on the agent positions 
        set_matrix(g,matrixW); % set g as a graph of matrixW
        connected=isconnected(g); %verify if g is a connected graph
        nd=nv(g);
        if connected            
            count=count+nrNodes;

            if strcmp(param.failureOp,'Random') % if random failures is considered
                
                position(randi([1 nd]),:)=[];  % remove from initial position matrix a random number (nd) of elements
                
            else
                if strcmp(param.failureOp,'BC') % for BC strategy
                    matrixW = initialize_matrixA(position,param,options); % initialized the adjacency matrix 
                    BC=[betweenness_wei(matrixW)  [1:nd].'];   % compute the node's BC
                    [ranking pos] = sort(BC(:,1),'descend');  % ordering nodes according to their BC ranking
                    position(pos(1:nrNodes),:)=[];  %remove higher BC nodes - the nunber of nodes is the same of vulnerable nodes
                end
            end
           
        end
 
    end
    
    if nv(g)==0 
        fnv=1; 
    else
        fnv=max(subgrafosSize(g))/nv(g);   % verifies the fraction of node into the giant component
    end
    
    % count/nN = the fraction of central nodes that need to be removed from the network to obtain a disconnected network. ]
    % count = number of nodes removed from the network
    % nN = number of nodes into the original (complete) network
    
    data=[count/nN fnv];

    free(g);
    
end
