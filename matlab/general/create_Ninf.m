% filename: create_Ninf.m
% Purpose:  creates a struct containing information (id and position) about the 2-hop neighborhood of each node 
% Input: 
% - position - x and y coordinates for each network agent
% - param - parametrization struct
% Output: 
% - data - cell of size (n,1), where n is the number of nodes into the network
% - each cell contains a cell of size (nN,3) where nN is the number of
% neighbors into that node + 1. 
% The first row encompasses: node id (col1) + node position (col 2) + list of node's neighbors
% For each of the node's neighbors there is a cell containing node id (col1) + node position (col 2) + list of node's neighbors
% Ex. consider node id 1 for data{1,1} - it contains is a cell of size {3,1}
% data{1,1}(1,1) = 1 - id
% data{1,1}(1,2) = [43.6 21.7] - node position
% data{1,1}(1,3) = [7 9] - node 1-hop neighbor
% data{1,1}(2,1) = 7 - id
% data{1,1}(2,2) = [37.6 42.7] - node position
% data{1,1}(2,3) = [1 8] - node 1-hop neighbor
% data{1,1}(2,1) = 9 - id
% data{1,1}(2,2) = [32.7 25.8] - node position
% data{1,1}(2,3) = [1 2] - node 1-hop neighbor


function [data] = create_Ninf(position,param)

matrix=initialize_matrixAdj(position,param); % create adjacency matrix

data=[];
for i=1:size(matrix,2)
    ddsN=[];
    listNi=kneighbors(matrix,i,1); % 1-hop neighbor
    ddsN={i position(i,:) listNi};
    for j=listNi
        ddsNj=[];
        listNj=kneighbors(matrix,j,1); % 2-hop neighbor
        for k=listNj
            if k~=i
                ddsNj=[ddsNj; k  position(k,:)];
            end
        end
        ddsN=[ddsN; {j} position(j,:) {ddsNj}];
    end
    data=[data; {ddsN}];
end