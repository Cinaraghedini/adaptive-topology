function [newC] = limitedVoronoiCell(position,param,options)  

distance = squareform(pdist(position,'euclidean'));

[m mw mwI] = initialize_matrix(distance,param,options);

voronoiPoints = pointVoronoiCell(positionNi);

