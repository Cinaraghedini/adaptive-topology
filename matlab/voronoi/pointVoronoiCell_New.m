function [newV newC]= pointVoronoiCell_New(position,param)

D = 10^5;

x=position(:,1);
y=position(:,2);

[vx,vy] = voronoi(x,y);

[V, C] = voronoin(position);

vx = round(vx*D)/D;
vy = round(vy*D)/D;
V = round(V*D)/D;

vxy=[vx(1,:)' vy(1,:)' vx(2,:)' vy(2,:)'];

vxyI=[vxy(:,1:2); vxy(:,3:4)];
u=unique(vxyI(:,1:2),'rows');

occurences=[];

for i=1:size(u,1)
    n=length(find(ismemberf(vxyI(:,1:2),u(i,1:2),'row', 'tol', param.tol)));
    occurences=[occurences;u(i,:) n];
end

occurences_u=find(occurences(:,3)==1);
if ~isempty(occurences_u)
    uniquePoints=occurences(occurences_u,1:2);
end

uniqueRef=[];
for i=1:size(uniquePoints,1)
    pos = find(ismemberf(vxy(:,1:2),uniquePoints(i,1:2),'row', 'tol',param.tol));
    if ~isempty(pos)
        uniqueRef=[uniqueRef; uniquePoints(i,1:2)  vxy(pos,3:4)];
    end
    pos = find(ismemberf(vxy(:,3:4),uniquePoints(i,1:2),'row', 'tol', param.tol));
    if ~isempty(pos)
        uniqueRef=[uniqueRef; uniquePoints(i,1:2)  vxy(pos,1:2)];
    end
end

refPoints=[];
newV=[];
newC=cell(size(position,1),1);

for i=1:size(C,1)
    point=V(C{i,1},:);
    newPoint=[];
    if ~isinf(point) %if there is no infinite points        
        [newPoint newV] = createCellPoints(point,newV,param);               
    else
        posInf=find(isinf(point(:,1))); %position of the inf point
        if posInf>1 && posInf<size(point,1) %inf is in the middle of the array points
            p1=point(posInf-1,:);
            p2=point(posInf+1,:);
            point=defineNewCellPoints(position(i,:),uniqueRef,point,p1,p2,posInf,'middle',param);   
            if ~isempty(point)
                [newPoint newV] = createCellPoints(point,newV,param);
            end
        else
            if posInf==1 %inf is in the first line
               p1=point(2,:); %reference point
               if (size(point,1)==2) %if there are only two points
                    point=defineNewCell2Points(position(i,:),uniqueRef,point,p1,posInf,'first');
                    if ~isempty(point)
                        [newPoint newV] = createCellPoints(point,newV,param);
                    end                    
                else %there are more than two point
                    p2=point(size(point,1),:); %reference point
                    point=defineNewCellPoints(position(i,:),uniqueRef,point,p1,p2,posInf,'first',param);
                    if ~isempty(point)
                        [newPoint newV] = createCellPoints(point,newV,param);
                    end
               end
            else % if inf point is in the last point position
                if posInf==size(point,1)
                    if (size(point,1)==2) %if there are only two points
                        p1=point(1,:); %first point is the reference point
                        point=defineNewCell2Points(position(i,:),uniqueRef,point,p1,posInf,'last');
                        if ~isempty(point)
                            [newPoint newV] = createCellPoints(point,newV,param);
                        end
                    else %if there are more than two points
                        p1=point(posInf-1,:);
                        p2=point(1,:);
                        point=defineNewCellPoints(position(i,:),uniqueRef,point,p1,p2,posInf,'last',param);
                        if ~isempty(point)
                            [newPoint newV] = createCellPoints(point,newV,param);
                        end                                                
                    end
                end                
            end            
        end
    end
    newC(i,1)={newPoint};
end