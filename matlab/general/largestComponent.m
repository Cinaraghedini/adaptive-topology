% filename: largestComponent.m
% Purpose:  find the largest connected component of a graph given by node
% positions and their range
% Input: 
% - position - x and y coordinates for each network agent
% - param - parametrization struct
% - options - graph options (set the parametrization for graph properties
% computation
% Output: 
% - positionS - matrix (n,2) where n is the number of  nodes, col 1 is the
% x position and col 2 is the y position into the plane


function positionS = largestComponent(position,param,options) 

positionS=position;
[m] = initialize_matrixAdj(position,param);
componentsG = components(m,options);
componentId=unique(componentsG);
k=[];
for i=1:size(componentId,1)
  k=[k;i length(find(ismember(componentsG,componentId(i,1))))];
end
largestComponent=find(k(:,2)==max(k(:,2)));
idS=k(largestComponent(1,1),1);
nodesS=find(~ismember(componentsG,idS));
positionS(nodesS,:)=[];