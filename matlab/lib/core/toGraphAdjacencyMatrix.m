%{
function toGraphAdjacencyMatrix: convert unsorted polygons vertexes to a graph
    represented as a adjacency matrix.

input: 
    unsortedPolygonsVertexes: unsorted edges of all possible polygons
        columns are similar to output polygons but last column represents 
        edges id.
    actualIndexes: GLOBAL vector mapping all actual indexes.

output: 
    graph adjacency matrix.
%}

function graph = toGraphAdjacencyMatrix(unsortedRobotsCenters)

    global actualIndexes tools;
    
    if(isempty(tools))
        tools = Tools();
    end

    [m, ~] = size(unsortedRobotsCenters);
    
    if(~isempty(unsortedRobotsCenters))
        actualIndexes = unique(unsortedRobotsCenters(:,3:4));
    else
        actualIndexes = [];
    end
    
    graph = zeros(length(actualIndexes));
 
    for i = 1 : m
        graph(tools.toRelativeIndex(unsortedRobotsCenters(i,3)), tools.toRelativeIndex(unsortedRobotsCenters(i,4))) = 1;    
    end

end