function [newC] = localLimitedVoronoiCell(position,param,options)  

distance = squareform(pdist(position,'euclidean'));

[m mw mwI] = initialize_matrix(distance,param,options);

for i=1:size(position,1)
 
        listNi=kneighbors(m,i,1);
        positionNi=[position(i,:)];
        for j=listNi
            positionNi=[positionNi;position(j,:)];
        end

        voronoiPoints = pointVoronoiCell(positionNi);
        
end