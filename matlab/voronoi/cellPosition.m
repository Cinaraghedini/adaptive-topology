% filename:cellPosition.m
% Purpose: verifies if a vertex is already into the matrix of vertices, if
% not, it is inserted into the mmatrix
% Input:
% - V - matrix of vertices 
% - newVertex - vertex to be evaluated 
% - param - parametrization struct
% Output:
% - V - matrix of vertices 
% - posVertex - the position of the vertex


function [V posVertex] = cellPosition(V,newVertex,param)

new=true;
if ~isempty(V)
    [~,indx]=ismemberf(newVertex,V,'row', 'tol', param.tol);
    if indx>0      
        posVertex=indx;
        new=false;
    end
end

if new
    V=[V;newVertex];
    posVertex=size(V,1);
end
       