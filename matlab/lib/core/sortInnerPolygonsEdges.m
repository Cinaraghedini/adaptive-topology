%{
function sortInnerPolygonsEdges: convert unsorted edges of robots - that belongs to a inner
    (blue) polygon composed by the bases of outer (red) intersections -
    into one or more closed polygons with vertexes correcty linked and ordered.
    
    The function selects the next child node based on the centroid created
    by all vertexes envolved (kept in a double linked list - the dlnode
    data structure).

input:
    unsortedRobotsCenters: unsorted edges of robots that belongs to a inner
        polygon composed by the bases of outer (red) intersections. This field
        is array having coordinate x of the robot at column 1, y at column
        2, the robot id at column 3, the other robot envolved at column 4,
        the cluster id at column 5 and edge id at column 6.
output:
    innerpolygons: global (see more on evalPolygons documentation).
        

@author Israel C. S. Rocha (israel.c.rocha@gmail.com)

%}

function sortInnerPolygonsEdges(graph)

    global robots tools globalGraph firstNode innerPolygons
    
    tools = Tools();
    polygonId = 0;    
    innerPolygons = [];

    globalGraph = graph;   
    
    while(true)
        nodesDegrees = sum(globalGraph);
        minDegrees = find(nodesDegrees == 2);

        if(isempty(minDegrees))
            break;
        end
        
        firstNode = selectFirstNodeByMinDegreeAndMaxCoordinate(minDegrees);

        performSort(firstNode);   

        polygonId = polygonId + 1;
        node = firstNode;
        firstRobotId = tools.toActualIndex(firstNode.Data);
        
        newInnerPolygon = [];
        while(~isempty(node))        
            robotId = tools.toActualIndex(node.Data);

            x = robots(robotId,1);
            y = robots(robotId,2);

            newInnerPolygon = [newInnerPolygon; x y robotId robotId, polygonId];

            node = node.Next;   
        end
        
        [m, ~] = size(newInnerPolygon);
        
        if(m <= 2)
            continue;
        end
        
        x = robots(firstRobotId, 1);
        y = robots(firstRobotId, 2);        
        innerPolygons = [innerPolygons; newInnerPolygon; x y firstRobotId firstRobotId polygonId]; %#ok<*AGROW>

    end
    
end

function performSort(fatherNode)

    global globalGraph firstNode
    
    father = fatherNode.Data;
    
    children = find(globalGraph(father, :) == 1);
    
    qtyOfChildren = length(children);    
    if(qtyOfChildren == 0)
        return;
    end    
        
    child = selectChildByCentroid(father, children);       
    globalGraph(father, child) = 0;
    globalGraph(child, father) = 0;
    
    if(child == firstNode.Data)
        return;
    end
    
    childNode = dlnode(child);
    childNode.insertAfter(fatherNode);

    performSort(childNode);
    
end

function child = selectChildByCentroid(father, children)

    if(length(children) == 1)
        child = children;
        return;
    end   

    global firstNode tools robots inner;
    
    centroid = findNodesCentroid();
    
    robotsIds = tools.toActualIndex(children);
    
    distances = pdist2(centroid, robots(robotsIds, 1:2));
    child = children(distances == min(distances));
    
    if(length(child) == 1)
        return;
    end
    
    fatherRobotId = tools.toActualIndex(father);
    
    envolvedInnerIntersections = inner((inner(:,3) == fatherRobotId | inner(:,4) == fatherRobotId), :);
    centroidOfEnvolvedInnerIntersections = mean(envolvedInnerIntersections(:,1:2));
    
    distances = pdist2(centroidOfEnvolvedInnerIntersections, robots(robotsIds, 1:2));
    child = child(distances == max(distances));
    
    if(length(child) == 1)
        return;
    end    
    
    robotFirst = robots(tools.toActualIndex(firstNode.Data), :);
    robotsChildren = robots(tools.toActualIndex(child), :);
    
    distances = pdist2(robotFirst(1:2), robotsChildren(:, 1:2));
    child = child(distances == min(distances));
    child = child(1);

end

function firstNode = selectFirstNodeByMinDegreeAndMaxCoordinate(minDegrees)

    global tools robots

    idsOfRobotsWithMinDegrees = tools.toActualIndex(minDegrees);        
    robotsXCoords = robots(idsOfRobotsWithMinDegrees, 1);
    idOfRobotWithMaxX = idsOfRobotsWithMinDegrees(robotsXCoords == max(robotsXCoords));

    if(length(idOfRobotWithMaxX) > 1)
        robotsXCoords = robots(idOfRobotWithMaxX, 2);
        idOfRobotWithMaxX = idsOfRobotsWithMinDegrees(robotsXCoords == max(robotsXCoords));
    end          

    firstNode = dlnode(tools.toRelativeIndex(idOfRobotWithMaxX(1)));

end

function point = findNodesCentroid()
    global firstNode tools robots;
    
    counter = 0;
    point = [0 0];
    
    node = firstNode;
    while(~isempty(node))        
        robotId = tools.toActualIndex(node.Data);
        point = point + robots(robotId, 1:2);
        counter = counter + 1;
        node = node.Next;   
    end
    
    point = point / counter;

end