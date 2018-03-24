% filename: attack_newtork.m
% Purpose:  removing N nodes from a network according to the strategy setting param.failureOp 
% Input: 
% - position - x and y coordinates for each network agent
% - param - parametrization struct
% - options - graph options (set the parametrization for graph properties
% computation
% Output: 
% - position - updated matrix of x and y coordinates for each network agent, i.e., without N nodes removed 
% - pos - list of nodes removed 


function [position pos]=attack_network(position,param,options)

nd=(size(position,1)); %number of nodes into the network

nrNodes=round(param.fractionIteration*nd); % number of nodes to be removed from the network

pos=[];

if nrNodes<nd || nd==1
    if strcmp(param.failureOp,'Random')
        count=0;
        while count<nrNodes
            pos=randi([1 nrNodes]); %generates a random number from a range
            position(pos,:)=[];  %remove from initial position matrix a random number of elements
            count=count+1;
        end
    else
        if strcmp(param.failureOp,'BC')
            matrixW = initialize_matrixA(position,param,options);
            BC=[betweenness_wei(matrixW)  [1:nd].'];   % compute the node's BC
            [ranking idx] = sort(BC(:,1),'descend');  % ordering nodes according to their BC ranking
            pos=idx(1:nrNodes);
            position(pos,:)=[];  %remove higher BC nodes - the nunber of nodes is the same of vulnerable nodes
        end
    end
end