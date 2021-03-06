% filename: defineNewCellPoints.m
% Purpose:  according to the vertices of the voronoi diagram replaces the
% Inf points in cells generated by voronoin with the best point, resulting
% in a voronoi cell without Inf points using as a reference the
% predecessors and sucessors vertices and the distance between them
% Input:
% - currentPosition: position of the generating point - node
% - uniqueRef: matrix with unique voronoi points 
% - verticesPositions: positions of vertives composing the voronoi cell for the generating point (node)
% - predecessorInf: predecessor of inf
% - successorInf: successor of inf
% - posInf: inf position into the vertices matrix
% - infID: 'middle', 'first' 'last' - label position of posInf
% - param - parametrization struct
% Output:
% - newVertices = returns the new vertices matrix without inf points
% Ex: 
% verticesPositions=[40.52221 21.829840;Inf Inf;43.8129 18.2939]
% newVertices=[40.52221 21.82984;53.09634 39.14236;58.52981 34.38280;43.81290 18.29392]

function newVertices=defineNewCellPoints(currentPosition,uniqueRef,verticesPositions,predecessorInf,successorInf,posInf,infID,param)   


newVertices=[];
candidate=[];
distP1_P2=[];
evaluatedPoints=[];
distP1_P2_evaluated=[];

posPredecessor=find(ismemberf(uniqueRef(:,3:4),predecessorInf,'row', 'tol', param.tol)); %finds if predecessor of inf is a unique point
posSuccessor=find(ismemberf(uniqueRef(:,3:4),successorInf,'rows', 'tol', param.tol)); %finds if successor of inf is a unique point
for ii=1:size(posPredecessor,1)
    refP1=uniqueRef(posPredecessor(ii),1:2); % set uniqueRef predecessor as the predecessor of inf  
    for k=1:size(posSuccessor,1)
        refP2=uniqueRef(posSuccessor(k),1:2); % set uniqueRef sucessor as the successor of inf 
        % creates matrix (aux) with the sequence of vertices according the inf
        % position
        if strcmp(infID,'middle') 
            aux=[verticesPositions(1:posInf-1,:);refP1;refP2;verticesPositions(posInf+1:size(verticesPositions,1),:)];
        else
            if strcmp(infID,'first')
                aux=[refP1; verticesPositions(2:size(verticesPositions,1),:);refP2];
            else
                if strcmp(infID,'last')
                    aux=[verticesPositions(1:posInf-1,:);refP1;refP2];
                end
            end
        end
        % verifies if the generating point (node) is inside of the polygono
        %  generated by the previous procedure. If yes, the sequence is
        %  candidate to replace the original
        if(inpolygon(currentPosition(1,1),currentPosition(1,2),aux(:,1),aux(:,2))) % True for points inside or on a polygonal region
            candidate=[candidate aux ];
            distP1_P2=[distP1_P2;pdist([refP1;refP2],'euclidean')]; 
        else % compute the distance between points evaluated and set them as evaluated 
            evaluatedPoints=[evaluatedPoints aux];
            d=abs(det([refP2-refP1;currentPosition-refP1]))/norm(refP2-refP1);
            distP1_P2_evaluated=[distP1_P2_evaluated;d];
        end
    end
end
if ~isempty(candidate) % if there are candidate vertices that creates a polygon by replacing inf points using the closest points
    posMax=find(distP1_P2==min(distP1_P2)); 
    posMax=(posMax*2)-1;
    newVertices=candidate(:,posMax:posMax+1);
else
    if ~isempty(evaluatedPoints) % search by the best point for replacing inf
        posMax=find(distP1_P2_evaluated==min(distP1_P2_evaluated));
        posMax=(posMax*2)-1;
        newVertices=evaluatedPoints(:,posMax:posMax+1);
    end
end    
