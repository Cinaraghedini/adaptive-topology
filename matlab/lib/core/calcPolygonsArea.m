%{
function calcPolygonsArea: calculates the outer polygons (red) and inner
    polygons (blue). see documentation of evalPolygons.

input:
    innerPolygon: global.
    outerPolygon: global.
output:
    area: total area of polygons

@author Israel C. S. Rocha (israel.c.rocha@gmail.com)
%}

function area = calcPolygonsArea()    
    area = removePolygonsInsidePolygons(performCalcInnerArea() + performCalcOuterArea());
end

%{
function performCalcInnerArea: calculates the inner
    polygons (blue). see documentation of evalPolygons.

input:
    innerPolygon: global.
    
output:
    area: total area of inner polygons
    innerPolygonsAreas: global matrix having:
        area of each inner polygon at column 1
        inner polygons identification at column 2

@author Israel C. S. Rocha (israel.c.rocha@gmail.com)
%}

function area = performCalcInnerArea()

    global innerPolygons innerPolygonsAreas    
    
    polygons = innerPolygons;    

    area = 0;
    if(isempty(polygons))
        return;
    end
    
    polygonIds = getUnique(polygons(:,5));
    
    [m, ~] = size(polygonIds);
           
    innerPolygonsAreas = nan(m, 2);
           
    for i = 1 : m        
        polygon = polygons(polygons(:,5) == polygonIds(i),:);                
        
        polygonArea = polyarea(polygon(:,1), polygon(:,2));
        area = area + polygonArea;
        
        innerPolygonsAreas(i, :) = [polygonArea, polygonIds(i)];        
    end
    
end

%{
function performCalcOuterArea: calculates the outer
    polygons (red). see documentation of evalPolygons.

input:
    outerPolygons: global.
    
output:
    area: total area of outer polygons
    
@author Israel C. S. Rocha (israel.c.rocha@gmail.com)
%}

function area = performCalcOuterArea()

    global outerPolygons    
    
    polygons = outerPolygons;    

    area = 0;
    if(isempty(polygons))
        return;
    end
    
    polygonIds = getUnique(polygons(:,5));
    
    [m, ~] = size(polygonIds);
           
    for i = 1 : m        
        polygon = polygons(polygons(:,5) == polygonIds(i),:);                
        
        polygonArea = polyarea(polygon(:,1), polygon(:,2));
        area = area + polygonArea;     
    end
    
end

%{
function removePolygonsInsidePolygons: evaluates each inner polygon to find
    and remove overlapping of area. If a inner polygon B is inside of
    another A, the area of polygon B is subtracted twice from overall area.

    innerPolygon: global.
    outerPolygon: global.
    overlapedInnerPolygons: global

output:
    area: true overall area of all polygons

@author Israel C. S. Rocha (israel.c.rocha@gmail.com)
%}

function area = removePolygonsInsidePolygons(area)

    global innerPolygons overlapedInnerPolygons innerPolygonsAreas
    
    overlapedInnerPolygons = [];
    
    [m, ~] = size(innerPolygonsAreas);
    if(isempty(innerPolygons) || m < 2)
        return;
    end
    
    sortedInnerPolygonsAreas = sortrows(innerPolygonsAreas, -1);
    
    polygonsIds = sortedInnerPolygonsAreas(:,2);
    
    for i = 1 : length(polygonsIds) - 1
        polygonAId = polygonsIds(i);
        polygonA = innerPolygons(innerPolygons(:,5) == polygonAId, :);
                
        for j = i + 1 : length(polygonsIds)
            polygonBId = polygonsIds(j);
            polygonB = innerPolygons(innerPolygons(:,5) == polygonBId, :);
            meanPoint = mean(polygonB(:,1:2));
            inside = inpolygon(meanPoint(1), meanPoint(2), polygonA(:,1), polygonA(:,2));
            if(~inside)
               continue;
            end
           
            area = area - 2*innerPolygonsAreas(innerPolygonsAreas(:,2) == polygonBId, 1);
            overlapedInnerPolygons = [overlapedInnerPolygons; innerPolygons(innerPolygons(:,5) == polygonBId, :)]; %#ok<AGROW>
            innerPolygons(innerPolygons(:,5) == polygonBId, :) = [];

        end
    end

end