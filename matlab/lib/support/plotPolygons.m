function plotPolygons(polygons, color, opaque)

     if(isempty(polygons))
        return;
     end
    
     if(nargin == 2 || ~exist('opaque', 'var'))
         opaque = false;
     end
    
    polygonIds = getUnique(polygons(:,5));
    
    [m, ~] = size(polygonIds);
       
    for i = 1 : m        
        polygon = polygons(polygons(:,5) == polygonIds(i), :);                
       
        h = fill(polygon(:,1), polygon(:,2), color);
        if(opaque)
            opacity = 0.7;
        else
            opacity = 0.3;
        end
        set(h,'facealpha', opacity);
    end

end