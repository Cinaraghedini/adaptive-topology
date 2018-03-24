%{
function evalIntersections: calculate intersections and classify each one as inner
or outer intersections.
input:
    robots: GLOBAL array having robots intances at lines
        x at column 1, y at 2
        and circular range (radius) at column 3.

output:
    inner: GLOBAL array having inner intersections
    outer: GLOBAL array having outer intersectios, that is, intersections at most
    external cluster countour.
    both arrays have structure like;
        column 1: x-coordinate of intersection
        column 2: y-coordinate of intersection
        column 3: robot index envolved
        column 4: another robot index envolved  
      
    Additionally the function adds a fourth column to robots global matrix
    indicating the cluster identification where the robot belongs to.

@author Israel C. S. Rocha (israel.c.rocha@gmail.com)

%}

function evalIntersections()

    global inner outer robots intersects;

    [m, ~] = size(robots);
    
    if(m == 1)
        robots(1, 4) = 1;
        return;
    end
    
    outer = [];
    inner = [];
    intersects = [];
        
    for i = 1 : m - 1
        for j = i + 1 : m
            
             robots(:, 4) = 1;
            clusterId = 0;

            if(isnan(i) || isnan(j) || i == j)                
                return;
            end

            % caculate intersection between robots and remove second
            % element when equal to first
            [x, y] = circcirc(robots(i,1), robots(i,2), robots(i,3), robots(j,1), robots(j,2), robots(j,3));

            if(x(1) == x(2) && y(1) == y(2))
                x(2) = [];
                y(2) = [];
            end        

            k = find(isnan(x));
            x(k) = [];
            y(k) = [];

            intersects = [intersects; [x', y', repmat(i, length(x),1), repmat(j, length(x),1)]];             %#ok<AGROW>
            
            
            % set cluster id to neighboors
            if(isnan(robots(i,4)) && isnan(robots(j,4)))
                clusterId = clusterId + 1;
                robots(i,4) = clusterId;
                robots(j,4) = clusterId;
            else if(isnan(robots(i,4)))
                    robots(i,4) = clusterId;
                else if(isnan(robots(j,4)))
                    robots(j,4) = clusterId;
                    end
                end
            end
        end
    end
    
    classifyIntersections();
    
    setRemainingClusterId();
            
end

function classifyIntersections()

    global inner outer intersects robots
   
        % calculate distances between intersections and all robots and subtract the robot's radius to
        % determine whether intersections is inner or outer, that is, if
        % the intersection is inside the most external border of the whole picture
        
        [m, ~] = size(intersects);
        
        if(m == 0)
            return;
        end
        
        distancesIntersecToCentersOfRobots = min(pdist2([robots(:,1), robots(:,2)], intersects(:, 1:2)) - repmat(robots(:,3), 1, m));
                
        i = find(lesserThan(distancesIntersecToCentersOfRobots, 0));
        o = setdiff(1 : m, i);
        
        inner = intersects(i, :);
        outer = intersects(o, :);       
end

function setRemainingClusterId()
    global clusterId robots
    
    clusterId = clusterId + 1;
    robots(isnan(robots(:,4)), 4) = clusterId;
    
end