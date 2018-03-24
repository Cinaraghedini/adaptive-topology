function [newV] = boundedVoronoiPoints(position,param)


% replacing infinitive vertices with real values

[V C] = pointVoronoiCell(position);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Defining cutting points for unbounded cells with one vertice inside the
% circle and another outside
% Defining edges lying inside the cell.


newV=[];

for j=1:size(C,1)
    center=position(j,:);
    xypoint=C{j,1};
    newPoints=[];
    for jj=1:size(xypoint,1)
        v1=V(xypoint(jj,1),:);
        v2=V(xypoint(jj,2),:);
        boundingSet = verify_boundingCircle(v1,v2,center,(param.range*0.8));
        outSidePoint=find(boundingSet==-1);
        if isempty(outSidePoint)
            newPoints=[newPoints; v1 v2] ;
        else
            if (length(outSidePoint)==1)
                insideP=find(boundingSet==1);
                posboundingInside=xypoint(jj, insideP);
                insideV=V(posboundingInside,:);
                [xl, yl]=compute_x(v1,v2,center(1,1),center(1,2),(param.range*0.8));
                posboundingOutside=xypoint(jj,outSidePoint);
                outsideV=V(posboundingOutside,:);
                [xy] = findIntersectionPoint([xl yl],outsideV);
                newPoints=[newPoints; insideV xy];
                 distancei_ii=pdist([outsideV; center],'euclidean');
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Creating new curve points for unbonded cells with both vertices outside
                % the range
                pointsCurve_Center= center+((outsideV-center)*(param.range*0.8))/distancei_ii;
                newPoints=[newPoints; xy pointsCurve_Center];
            else
                distanceP_center=pdist([v1; center],'euclidean');
                newPointV1 = [center+((v1-center)*(param.range*0.8))/distanceP_center];
                distanceP_center=pdist([v2; center],'euclidean');
                newPointV2 = [center+((v2-center)*(param.range*0.8))/distanceP_center];
                newPoints=[newPoints; newPointV1 newPointV2];
            end
        end
        
    end
    
    
    auxEdgePoints=([newPoints(:,1:2);newPoints(:,3:4)]);
    edgePoints=[];
    for ii = 1:length(auxEdgePoints)
        indx=length(find(ismember(auxEdgePoints(:,1:2),auxEdgePoints(ii,1:2),'rows')));
        if indx==1
            edgePoints=[edgePoints; auxEdgePoints(ii,:)];
        end
    end
    
    if ~isempty(edgePoints)
         newPoints=[ newPoints; edgePoints(1,:) edgePoints(2,:)];
    end
        
    newV=[newV; newPoints];

end
