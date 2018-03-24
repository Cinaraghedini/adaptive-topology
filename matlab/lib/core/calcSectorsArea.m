%{
function calcSectorsArea: calculates the remaining (yellow) area not covered by
    outer or inner polygons. These sectors are composed by circles and its
    outer intersection points not present in its outer polygons.
input:
    robots: GLOBAL array having robots intances at lines
        x at column 1, y at 2
        and circular range (radius) at column 3.
    outer: GLOBAL (see evalIntersection documentation)
    inner: GLOBAL (see evalIntersection documentation)
    outerPolygons: GLOBAL (see evalPolygons documentation)

output:
    area: total area of all sectors in the field.

@author Israel C. S. Rocha (israel.c.rocha@gmail.com)

%}

function area = calcSectorsArea()

    global sectors sectorId
    
    sectorId = 0;
    sectors = [];
                       
    area =  getNotClusteredSectorsArea + getClusteredSectorsArea();
   
end

function area = getNotClusteredSectorsArea()

    global robots

    area = 0;    
    
    idsOfRobotWithoutIntersections = getIdsOfRobotWithoutIntersections();
    
    for i = 1 : length(idsOfRobotWithoutIntersections)
        robot = robots(idsOfRobotWithoutIntersections(i),:);
        R = robot(3);
        center = [robot(1), robot(2)];
        buildSector(center, R, 0, 360);
        area = area + getSectorArea(R, 0, 360);
    end;

end

function area = getClusteredSectorsArea()

    global outer robots
    
    area = 0;
    
    clustersIds = getUnique(robots(:, 4));
    if(~isempty(outer))
        outerIds = getUnique(outer(:,3:4));
    else
        outerIds = [];
    end
    
    for i = 1 : length(clustersIds)
    
        robotsIds = find(robots(:,4) == clustersIds(i));
        robotsIds = intersect(robotsIds, outerIds);
        
        for j = 1 : length(robotsIds)
            robotId = robotsIds(j);            
            robotOuters = outer((outer(:,3) == robotId | outer(:,4) == robotId), :);
            area = area + getSectorsArea(robotId, robotsIds, robotOuters);
        end
    
    end
    
end

function area = getSectorsArea(robotId, robotsIds, robotsOuters)

    global robots

    center = [robots(robotId,1) robots(robotId,2)];
    R = robots(robotId, 3);
    
    [m, ~] = size(robotsOuters);    
    if(m == 0 || m == 1)
        buildSector(center, R, 0, 360);
        area = getSectorArea(R, 0, 360);
        return;
    end
    
    area = 0;
    
    alfas = sort(getPolarCoordinatesRelativeToRobotCenter(robotsOuters, center));
    
    for i = 1 : length(alfas)
       
        j = i + 1;
        if(j > length(alfas))
            j = mod(j, length(alfas));
        end
        
        alfa = alfas(i);
        beta = alfas(j);
        
        if(beta < alfa)
            beta = 360 + beta;
        end
        
        if(isInsideSector(alfa, beta, center, R, robotsIds))
            continue;
        end
        
        buildSector(center, R, alfa, beta);        
        area = area + getSectorArea(R, alfa, beta);
        
        
    end
    
end

function alfas = getPolarCoordinatesRelativeToRobotCenter(points, center)

    [mr, ~] = size(points);
    alfas = NaN * ones(1, mr);
    
    for i = 1 : mr
        vertex = points(i,:);        
                   
        x = vertex(1);        
        y = vertex(2);        
       
        relativeX = x - center(1);
        relativeY = y - center(2);

        alfas(i) = getAngle(relativeX, relativeY);                
    end
    
end

function inside = isInsideSector(fromAngle, toAngle, actualRobotCenter, actualRobotRadius, robotsInClusterIds)

    global robots
    
    delta = toAngle - fromAngle;

    epsilon = delta / 100;
    d = (fromAngle + epsilon) * pi / 180;

    xd = actualRobotRadius * cos(d) + actualRobotCenter(1);
    yd = actualRobotRadius * sin(d) + actualRobotCenter(2);
    
    distances = pdist2([xd, yd], [robots(robotsInClusterIds,1), robots(robotsInClusterIds,2)]) - robots(robotsInClusterIds,3)';
    distances = min(distances);
    
    inside = lesserThan(min(distances), 0);
end

function buildSector(center, R, alfa, beta)

    global sectors sectorId

    [X, Y] = getCircleSector(center, R, alfa, beta);
       
    sectorId = sectorId + 1;
    Z = sectorId * ones(length(X),1);

    sectors = [sectors; X', Y', Z];

end

function area = getSectorArea(R, alfa, beta)
     angle = beta - alfa;
     area = pi * power(R, 2) * angle / 360;
end


function angle = getAngle(x, y)
    angle = atan2d(y, x);  
    
    if(angle < 0)
        angle = 360 + angle;
    end
end

function idsOfRobotWithoutIntersections = getIdsOfRobotWithoutIntersections()

    global robots outer inner
    
    [m, ~] = size(robots);

    if(~isempty(outer))
        robotsIdsThatHasOuterIds = getUnique(outer(:, 3:4));        
    else
        robotsIdsThatHasOuterIds = [];
    end
    
    idsOfRobotWithoutIntersections = setdiff(1 : m, robotsIdsThatHasOuterIds);
    
    
    if(~isempty(inner))
        robotsIdsThatHasOuterIds = getUnique(inner(:, 3:4));        
    else
        robotsIdsThatHasOuterIds = [];
    end
    
    idsOfRobotWithoutIntersections = setdiff(idsOfRobotWithoutIntersections, robotsIdsThatHasOuterIds);

end