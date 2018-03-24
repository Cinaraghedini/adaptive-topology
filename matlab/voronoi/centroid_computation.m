
function [centroidN centroid]=centroid_computation(x)

% matrix with the position x,y  
% centroid position x,y of the centroidal position 

centroid=[];
centroidN=[];
%Extract topological information of the diagram

[v,c] = voronoin([x(:,1) x(:,2)]);

for i=1:length(c)    
    Xi = v(c{i},:);
    centroid(i,1)=-1;
    centroid(i,2)=-1;
    if ~isinf(Xi)
        [geom,iner,cpmo]=polygeom(Xi(:,1),Xi(:,2));
        centroid(i,1)=geom(2);
        centroid(i,2)=geom(3);        
        centroidN=[geom(2) geom(3)];        
    end    
end


