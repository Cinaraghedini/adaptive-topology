function [centroids]=define_centroid(position,param)


centroids=position;

close all;

[newCells]=boundedVoronoi(position,param); %bounded voronoi cells

% figure(1)
% 
% vx=[];
% vy=[];
% 
% for ii=1:size(newCells,1)
%     cellPoints=newCells{ii,1};
%     for jj=1:size(cellPoints,1)
%         vx=[vx [cellPoints(jj,1);cellPoints(jj,3)]];
%         vy=[vy [cellPoints(jj,2);cellPoints(jj,4)]];
%     end
% end
% 
% plot(position(:,1),position(:,2),'r+',vx,vy,'b-');
% 
% 
% figure(2)
% 
% x=position(:,1);
% y=position(:,2);
% [vx,vy] = voronoi(x,y);
% 
% plot(x,y,'r+',vx,vy,'b-');


for i=1:size(newCells,1)
    
    [XY]=createVoronoiPoints(newCells{i,1},position,param);
    
    if ~isempty(XY)
        [geom,iner,cpmo]=polygeom(XY(:,1),XY(:,2));
        centroids(i,1)=geom(2);
        centroids(i,2)=geom(3);    
    end
    
end

%disp('leave defined_centroid')