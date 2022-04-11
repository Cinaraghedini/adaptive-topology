% filename: pointVoronoiCell.m
% Purpose:  generates a voronoi diagram from a set of generating points
% without infinite vertices
% Input:
% - position - matrix with x and y coordinates for each generating point
% for a voronoi diagram
% - param - parametrization struct
% Output: newCells contain one entry for each generating point. newVertices
% contains the total voronoi diagram vertices.
% newVertices = new vertices of a voronoi diagram after preprocessing de inf values replacing them for non infinite values
% newCells = cell contaning the sequence of vertices in order to generate a bounded voronoi cell
% - newVertices= [40.522210 21.829840;53.096340 39.142360;58.529810 34.382800;43.812900 18.293920]
% - newCells : [1 2;2 3;3 4;4 1]

function [newVertices, newCells]= pointVoronoiCell(position,param)

D = 10^5; % for rounding the voronoi edges 

x=position(:,1);
y=position(:,2);

[vx,vy] = voronoi(x,y); % returns bounded cells of the Voronoi diagram for nodes positioned at coordinates x and y

vx = round(vx*D)/D; %rounds voronoi cell vertices
vy = round(vy*D)/D; %rounds voronoi cell vertices

vxy=[vx(1,:)' vy(1,:)' vx(2,:)' vy(2,:)']; % transposes vertices of the Voronoi cells vertices

vxyI=[vxy(:,1:2); vxy(:,3:4)]; %creates a matrix with the Voronoi cells vertices at x and y coordinates

u=unique(vxyI(:,1:2),'rows'); % finds the unique vertices into the Voronoi diagram

[V, C] = voronoin(position); % returns Voronoi vertices V and the Voronoi cells C 

V = round(V*D)/D; %rounds voronoi cell vertices

occurences=[];

% computese the number of occurences of each vertice
for i=1:size(u,1)
    n=length(find(ismemberf(vxyI(:,1:2),u(i,1:2),'row', 'tol', param.tol)));
    occurences=[occurences;u(i,:) n];
end

% finds vertices that are unique, i.e., it is not connecting two vertices and
% thus, the cell is unbouded 

occurences_u=find(occurences(:,3)==1);
if ~isempty(occurences_u)
    uniquePoints=occurences(occurences_u,1:2);
end


% for each unique vertices finds the vertice it is connected to
uniqueRef=[]; 
for i=1:size(uniquePoints,1) 
    pos = find(ismemberf(vxy(:,1:2),uniquePoints(i,1:2),'row', 'tol',param.tol)); 
    if ~isempty(pos)
        uniqueRef=[uniqueRef; uniquePoints(i,1:2)  vxy(pos,3:4)];
    end
    pos = find(ismemberf(vxy(:,3:4),uniquePoints(i,1:2),'row', 'tol', param.tol));
    if ~isempty(pos)
        uniqueRef=[uniqueRef; uniquePoints(i,1:2)  vxy(pos,1:2)];
    end
end

% vx,vy generated by voronoi function returns the voronoi diagram vertices without the information about the node belonging to each cell. 
% On the other hand, voronoin generates each node cell information but with infinite vertices. Sometimes the inf is representing 
% more than one vertices. Next code block aimis at replacing each of the inf vertices in V to a vertice in vx,vy.    

newVertices=[];
newCells=cell(size(position,1),1);


for i=1:size(C,1) % For each voronoi cell
    verticesPosition=V(C{i,1},:); % vertices composing the i voronoi cell
    newPoint=[];
    if ~isinf(verticesPosition) %if there is no infinite points        
        [newPoint newVertices] = createCellPoints(verticesPosition,newVertices,param); % creates/updates cell for node i              
    else % if there is points to inf
        posInf=find(isinf(verticesPosition(:,1))); %inf position into the vertices matrix
        if posInf>1 && posInf<size(verticesPosition,1) % inf is in the middle of the vertices matrix
            predecessorInf=verticesPosition(posInf-1,:); % predecessor of inf
            successorInf=verticesPosition(posInf+1,:); % successor of inf
            verticesPosition=defineNewCellPoints(position(i,:),uniqueRef,verticesPosition,predecessorInf,successorInf,posInf,'middle',param); % replace inf vertices with non inf potin  
            if ~isempty(verticesPosition) % if it was possible to find a polygon with the available information set it as newPoint 
                [newPoint newVertices] = createCellPoints(verticesPosition,newVertices,param);
            end
        else 
            if posInf==1 %inf is in the first line
               predecessorInf=verticesPosition(2,:); %reference point
               if (size(verticesPosition,1)==2) %if there are only two vertices and one inf
                    verticesPosition=defineNewCell2Points(position(i,:),uniqueRef,verticesPosition,predecessorInf,posInf,'first');
                    
                 %   defineNewCell2Points(position,uniqueRef,verticesPosition,nonInfPosition,posInf,infID)
                    
                    if ~isempty(verticesPosition) % if it was possible to find a polygon with the available information set it as newPoint 
                        [newPoint newVertices] = createCellPoints(verticesPosition,newVertices,param);
                    end                    
                else %there are more than two point
                    successorInf=verticesPosition(size(verticesPosition,1),:); %reference point
                    verticesPosition=defineNewCellPoints(position(i,:),uniqueRef,verticesPosition,predecessorInf,successorInf,posInf,'first',param);
                    if ~isempty(verticesPosition)
                        [newPoint newVertices] = createCellPoints(verticesPosition,newVertices,param);
                    end
               end
            else % if inf point is in the last point position
                if posInf==size(verticesPosition,1)
                    if (size(verticesPosition,1)==2) %if there are only two  vertices and one inf
                        predecessorInf=verticesPosition(1,:); %first point is the reference point
                        verticesPosition=defineNewCell2Points(position(i,:),uniqueRef,verticesPosition,predecessorInf,posInf,'last');
                        if ~isempty(verticesPosition)
                            [newPoint newVertices] = createCellPoints(verticesPosition,newVertices,param);
                        end
                    else %if there are more than two points
                        predecessorInf=verticesPosition(posInf-1,:);
                        successorInf=verticesPosition(1,:);
                        verticesPosition=defineNewCellPoints(position(i,:),uniqueRef,verticesPosition,predecessorInf,successorInf,posInf,'last',param);
                        if ~isempty(verticesPosition)
                            [newPoint newVertices] = createCellPoints(verticesPosition,newVertices,param);
                        end                                                
                    end
                end                
            end            
        end
    end
    newCells(i,1)={newPoint};
end