%{
function evalPolygons: finds and delimits each outer and inner polygon
input:
    inner: GLOBAL output of function evalIntersections
    outer: GLOBAL output of function evalIntersections
    robots: GLOBAL output of function evalIntersections

output:
    outerPolygons: (red) triangles having outer points e robots centers as
    vertexes.
    innerPolygons: (blue) polygons formed by bases of outerPolygons 
    
    Both are matrixes having
        x-coordinate of polygon vertex  at column 1
        y-coordinate of polygon vertex  at column 2
        robot id envolved at column 3
        another robot id envolve at column 4
        polygon's id at column 5.

@author Israel C. S. Rocha (israel.c.rocha@gmail.com)

%}

function evalPolygons()
    
    evalOuter();
    evalInner();

end

%{
function evalOuter: evaluates (red) outer triangles composed by outer
    intersections and the geometric centers of two envolved robots.

input: 
    robots: GLOBAL
    outer: GLOBAL

output:
    outerPolygons: array having vertex x at column 1, vertex y at column 2, 
        robot envolved id at column 3, second robot envolved at column 4
        and polygon id at column 5.
%}

function evalOuter()
    
    global robots outer outerPolygons
    
    polygonOrder = 0;
    outerPolygons = [];
    triangleBases = [];
    
    [m, ~] = size(robots);
    if(m < 2)
        return;
    end
    
    [m, ~] = size(outer);
        
    if(m == 0)
        return;
    end       
    
    for i = 1 : m
                
        robotOrig = outer(i, 3);
        robotDest = outer(i, 4);
                
        polygonOrder = polygonOrder + 1;
                
        triangleBase = [
                        getPoint(robots(robotOrig, 1), robots(robotOrig,2), robotOrig, robotDest, polygonOrder)
                        getPoint(robots(robotDest, 1), robots(robotDest,2), robotDest, robotOrig, polygonOrder)
                        ];
        
        newPolygon = [
                getPoint(outer(i,1), outer(i,2), robotOrig, robotDest, polygonOrder)
                triangleBase
                getPoint(outer(i,1), outer(i,2), robotOrig, robotDest, polygonOrder)
                ];
            
        triangleBases = [triangleBases; triangleBase]; %#ok<AGROW>
        
       if(isCollateralVertexes(newPolygon))
           continue;
       end
       
        outerPolygons = [outerPolygons; newPolygon]; %#ok<AGROW>
        
    end
    
end

%{
function innerOuter: evaluates (blue) inner polygons from bases of outer polygons and
    innner intersections.

input: 
    robots: GLOBAL
    inner: GLOBAL
    outerPolygons: GLOBAL

output:
    innerPolygons: array having vertex x at column 1, vertex y at column 2, 
        robot envolved id at column 3, second robot envolved at column 4
        and polygon id at column 5.
%}

function evalInner()

   global outerPolygons polygonOrder
   
   if(isempty(outerPolygons))
       return;
   end
      
   polygonOrder = 0;
  
   polygonIds = getUnique(outerPolygons(:,5));
         
   unsortedRobotsEdges = [];
   for i = 1 : length(polygonIds)

       k = find(outerPolygons(:, 5) == polygonIds(i));
       if(isempty(k))
           continue;
       end
       edgesOfTriangle = outerPolygons(k, :);

       if(isCollateralVertexes(edgesOfTriangle))
           continue;
       end

       vertexA = edgesOfTriangle(2,:);
       vertexB = edgesOfTriangle(3,:);
       triangleBase = edgesOfTriangle(2:3, :);       
       idOfRobotA = vertexA(3);
       idOfRobotB = vertexB(3);
             
       if(isNotInSameCluster(idOfRobotA, idOfRobotB))
           continue;
       end
              
       if(isEdgeInTwoOuterPolygons(idOfRobotA, idOfRobotB))
           continue;
       end
       
       unsortedRobotsEdges = [unsortedRobotsEdges; triangleBase]; %#ok<AGROW>
       
   end

   sortInnerPolygonsEdges(toGraphAdjacencyMatrix(unsortedRobotsEdges));
   
   
end

%{
function isCollateralVertexes: whether triangle's vertexes are collateral,
    that is, zero area triangle.

input: 
    precision: GLOBAL decimal precision.
    triangle: vertexes of the triangle.

output:
    boolean indicanting wheter is collateral.
%}

function collateral = isCollateralVertexes(triangle)

    D = ones(3);
    D(1:3, 1:2) = triangle(1:3, 1:2);
    
    collateral = isEqualTo(det(D), 0);
end

function ret = isNotInSameCluster(idOfRobotA, ifOfRobotB)

    global robots
    
    ret = robots(idOfRobotA, 4) ~= robots(ifOfRobotB, 4);

end

function ret = isEdgeInTwoOuterPolygons(idOfRobotA, idOfRobotB)

    global outerPolygons

     k = find((outerPolygons(:,3) == idOfRobotA & outerPolygons(:,4) == idOfRobotB));
     polygonsIds = getUnique(outerPolygons(k, 5)); %#ok<*FNDSB>
     
     ret = length(polygonsIds) ~= 1;
       
end