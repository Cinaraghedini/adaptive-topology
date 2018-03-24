% filename:createCellPoints.m
% Purpose: returns updated newV for the vertices (vertices) and a matrix
% linking vertices to define a bounded cell.
% Input:
% - vertices - vertex(ices) to be evaluated 
% - newV - vertices matrix of vertices 
% - param - parametrization struct
% Output:
% - inkedVertices - matrix containing the position of each vertex into the
% vertices matrix (newV) chaining the vertices sequence:
% linkedVertices= [1 2;2 3;3 4;4 1]
% - newV - updated vertices points [40.52221 21.82984;53.09634 39.14236; 58.52981 34.38280; 43.8129 18.29392];


function [linkedVertices newV] = createCellPoints(vertices,newV,param)

linkedVertices=[];
posAux=vertices(1,:); % first vertex 
[newV posVertex_A] = cellPosition(newV,posAux,param); % returns the vertex posAux position into the vertice matrix. If it is not into the matrix it is inserted 
for l=2:size(vertices,1) % for each of the remaining vertices 
    [newV posVertex_B] = cellPosition(newV,vertices(l,:),param); 
    linkedVertices=[linkedVertices;[posVertex_A posVertex_B]];
    posVertex_A=posVertex_B;
end
linkedVertices=[linkedVertices;[posVertex_A linkedVertices(1,1)]]; % connectining first and last vertices