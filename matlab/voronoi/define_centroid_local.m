% filename: define_centroid_local.m
% Purpose:  compute the centroid of each the voronoi cell using only local information
% Input:
% - position - matrix with x and y coordinates for each network agent
% - param - parametrization struct
% - options - graph options (set the parametrization for graph properties
% computation
% Output:
% - centroids - matrix with x and y coordinates for each network agent cell centroid

function [centroids]=define_centroid_local(position,options,param)

centroids=position;

distance = squareform(pdist(position,'euclidean')); % computes the distance between nodes 

[m mw mwI] = initialize_matrix(distance,param,options); % initializes the network matrices

close all;

figure(1)

x=position(:,1);
y=position(:,2);

figure(1)

[vx,vy] = voronoi(x,y);

plot(x,y,'b*',vx,vy,'k.-');

hold all;

for j=1:size(position,1)
    disp(j)
    listNi=kneighbors(m,j,1); % finds the neighbors of node j
    
    if length(listNi)>1 %if node j has any neighbors
        positionNi=position(j,:); % initialize a matrix that will contain the position of node j itself, as well as its neighbors
        for jj=listNi
            positionNi=[positionNi;position(jj,:)];
        end
        [newCells]=boundedVoronoi(positionNi,param); % generates a bounded voronoi diagram for node j local network  
        [XY]=createVoronoiPoints(newCells{1,1},param); % generates the polygon for node j voronoi cell 
        if ~isempty(XY)
            [geom,iner,cpmo]=polygeom(XY(:,1),XY(:,2)); % returns area, X centroid, Y centroid and perimeter for the planar polygon specified by vertices in vectors X and Y.
            centroids(j,1)=geom(2); % set x coordinate for node j cell centroid
            centroids(j,2)=geom(3); % set y coordinate for node j cell centroid
        end
        gcf1=plotTestVoronoi(positionNi,centroids,newCells,j);
       
       vx1=[];
       vy1=[];
       
       cellPoints=newCells{1,1};
       for jj=1:size(cellPoints,1)
           vx1=[vx1 [cellPoints(jj,1);cellPoints(jj,3)]];
           vy1=[vy1 [cellPoints(jj,2);cellPoints(jj,4)]];
       end
       
       plot(vx1,vy1,'LineStyle',':','LineWidth',1,'Color',[0.800000011920929 0 0.800000011920929]);
       
       hold all;
       
       plot(centroids(1,1),centroids(1,2),'g*');

       hold all;
       
    end
end