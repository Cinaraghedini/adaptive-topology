% filename: createVoronoiPoints.m
% Purpose:  converting the edges in the format generated by the
% preprocessing of voronoi diagram to a matrix with sequencially ordered
% vertices. The vertices are processed in "clock-wise". For each vertex, the vertex on its "right" is searched. Notice that each vertex 
% has exactly two vertices, thus being at two rows (because a preprocessed voronoi cells is bounded)   
% Input:
% - voronoiPoints: cell containing the vertices of a voronoi cell. Each
% edge is a matrix row with its two vertices.
% newC{1,1}(1,1) x coordinate of the first vertex, 
% newC{1,1}(1,2) y coordinate of the first vertex,  
% newC{1,1}(1,3) x coordinate of the second vertex
% newC{1,1}(1,4) y coordinate of the secibd vertex.  
% - param - parametrization struct
% Output:
% - XY: matrix containing the x and y coordinates of the each voronoi
% vertex sequencialy ordered.


function [XY]=createVoronoiPoints(voronoiPoints,param)

XY=[voronoiPoints(1,1:2); voronoiPoints(1,3:4)]; % first row is the  reference for starting the vertices list

ref=voronoiPoints(1,3:4);  % defines reference points as the one at the right (clock-wise) of the edge, using first row 

refPos=1;

nrPoints=size(unique([voronoiPoints(:,1:2); voronoiPoints(:,3:4)],'rows'),1); % number of unique vertices

while(size(XY,1)<nrPoints) % while there are vertices to be processed
   
    found=false;
    
    pos=find(ismemberf(voronoiPoints(:,1:2),ref,'rows', 'tol', param.tol)); % finds vertices into the original cell taking left positions as reference
    
    % if the reference vertex is at the left position in row
    if ~isempty(pos)  % if node is in more than one row
        if length(pos)>1 % removes the refPos index from pos 
            pos=setdiff(pos,refPos); 
        end
        if pos~=refPos % if vertex position and refPos are not the same vertex the right position is added to the vertices lists  
            XY=[XY;voronoiPoints(pos,3:4)];
            ref=voronoiPoints(pos,3:4); % the new vertex is now the reference for finding its compouding edge vertex
            refPos=pos; 
            found=true; % the reference vertex found its next point
        end
    end
    
    % if no vertex is found at the left side of the cell row, the vertex at the right side is now used for searching and the same process is applied 
    
    if ~found 
        pos=find(ismemberf(voronoiPoints(:,3:4),ref,'rows', 'tol', param.tol));
        if ~isempty(pos)
            if length(pos)>1
                pos=setdiff(pos,refPos);
            end
            if pos~=refPos % if vertex position and refPos are not the same vertex the left position is added to the vertices lists  
                XY=[XY;voronoiPoints(pos,1:2)];  % the new vertex is now the reference for finding its compouding edge vertex
                ref=voronoiPoints(pos,1:2); 
                refPos=pos;
                found=true;
            end
        end
    end    
end