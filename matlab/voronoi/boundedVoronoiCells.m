function [newC] = boundedVoronoiCells(position,param,options)  

distance = squareform(pdist(position,'euclidean'));

[m mw mwI] = initialize_matrix(distance,param,options);

[V C] = pointVoronoiCell(position);

vx=[];
vy=[];

x=position(:,1);
y=position(:,2);
 
figure(1)

for j=1:size(C,1)
    xypoint=C{j,1};
    for jj=1:size(xypoint,1)
        vx=[vx [V(xypoint(jj,1),1);V(xypoint(jj,2),1)]];
        vy=[vy [V(xypoint(jj,1),2);V(xypoint(jj,2),2)]];
    end
end

plot(x,y,'r+',vx,vy,'b-');

xlimit=[min(x)-20 max(x)+20];
ylimit=([min(y)-20 max(y)+20]);

xlim(xlimit);
ylim(ylimit);

saveas(gcf,[param.pathR,'preProcessed','.fig'],'fig');
saveas(gcf,[param.pathR,'PreProcessed','.png'],'png');


figure(2)

[vx,vy] = voronoi(x,y);

plot(x,y,'r+',vx,vy,'b-');

xlimit=[min(x)-20 max(x)+20];
ylimit=([min(y)-20 max(y)+20]);

xlim(xlimit);
ylim(ylimit);

saveas(gcf,[param.pathR,'Original','.fig'],'fig');
saveas(gcf,[param.pathR,'Original','.png'],'png');


new_V=[];
new_V_ext=[];
for i=1:size(position,1)
     pos_i=position(i,:);
     V_i=[];
     for j=1:size(C{i,1},1)
         V_i = [V_i; V(C{i,1}(j,1),:); V(C{i,1}(j,2),:)]; 
     end
     V_i= unique(V_i,'rows');
     for ii=1:size(V_i,1)
         pos_ii=V_i(ii,:);
         distancei_ii=pdist([pos_ii; pos_i],'euclidean');
         if distancei_ii > (param.range*0.8)
             [~,posVertex] = cellPosition(V,pos_ii,param);             
             new_V=[new_V; pos_i+((pos_ii-pos_i)*(param.range*0.8))/distancei_ii];
             new_V_ext=[new_V_ext; pos_i pos_ii pos_i+((pos_ii-pos_i)*(param.range*0.8))/distancei_ii];
              
         end
     end     
end

[xc,yc,r,a] = circfit(new_V(:,1),new_V(:,2));

[xCircle, yCircle ]=GeneratePointsCircle(r,xc,yc);


newC=cell(size(C,1),1);
newV=[];



cutPoint=[];

partialOutside=[];
completeOutside=[];

for j=1:size(C,1)
    
    xypoint=C{j,1};    
    newPoint=[];
    for jj=1:size(xypoint,1)
        v1=V(xypoint(jj,1),:);
        v2=V(xypoint(jj,2),:);
        boundingSet = verify_boundingCircle(v1,v2,[xc,yc],r);
        outSidePoint=find(boundingSet==-1);
        if isempty(outSidePoint)
            [newV posVertex_A] = cellPosition(newV,v1,param);
            [newV posVertex_B] = cellPosition(newV,v2,param);
            newPoint=[newPoint; posVertex_A posVertex_B];
        else
            
            if ~isempty(outSidePoint)
                
                if (length(outSidePoint)==1)
                    
                    inSidePoint=find(boundingSet==1);
                    
                    posbounding=xypoint(jj,inSidePoint);
                    
                    insideV=V(posbounding,:);
                    
                    [xl, yl]=compute_x(v1,v2,xc,yc,r);
                    
                    [xy] = findIntersectionPoint([xl yl],insideV);
                    
                    cutPoint=[cutPoint;xy];
                    
                    pointOutside=xypoint(jj,outSidePoint);
                    
                    outsideV=V(pointOutside,:);
                    
                    [newV posVertex_A] = cellPosition(newV,insideV,param);
                    [newV posVertex_B] = cellPosition(newV,xy,param);
                    
                    newPoint=[newPoint; posVertex_A posVertex_B];
                    
                                  
                   % partialOutside=[partialOutside; j jj v1 v2 posbounding insideV pointOutside outsideV xy];
                    
                                       
                else
                    
                    [xl, yl]=compute_x(v1,v2,xc,yc,r);
                    
                    [xy] = findIntersectionPoint([xl yl],v1);
                    
                    
                    
                    
    %                [xunit,yunit] = createSegmCircle(v1,v2,center);
                    
                    
                    %completeOutside=[completeOutside; j jj v1 v2];
                  
                end
            end
        end
    end
    
    uniquePoints=unique([newPoint(:,1); newPoint(:,2)]);
    uniquePoint=[];
    for i=1:size(uniquePoints,1)
        cellPoints=find([newPoint(:,1); newPoint(:,2)]==uniquePoints(i,1));
        if size(cellPoints,1)==1
           uniquePoint=[uniquePoint uniquePoints(i,1)];
        end
    end
        
    newC(j,1)={[newPoint; uniquePoint]};
end


figure(3)

plot(x,y,'r+',vx,vy,'b-');
hold all;
plot(new_V(:,1),new_V(:,2),'g*');
hold all;
plot(xCircle,yCircle,'g');
hold all;
plot(xc,yc,'mo');

xlimit=[min(x)-20 max(x)+20];
ylimit=([min(y)-20 max(y)+20]);

xlim(xlimit);
ylim(ylimit);

saveas(gcf,[param.pathR,'preProcessed_2','.fig'],'fig');
saveas(gcf,[param.pathR,'PreProcessed_2','.png'],'png');


figure(4)

plot(x,y,'r+',vx,vy,'b-');
hold all;
plot(new_V(:,1),new_V(:,2),'g*');
hold all;
plot(xCircle,yCircle,'k+');
hold all;
plot(xc,yc,'mo');
hold all;
plot(cutPoint(:,1),cutPoint(:,2),'csquare');
hold all;


xlimit=[min(x)-20 max(x)+20];
ylimit=([min(y)-20 max(y)+20]);

xlim(xlimit);
ylim(ylimit);


saveas(gcf,[param.pathR,'preProcessed_3','.fig'],'fig');
saveas(gcf,[param.pathR,'PreProcessed_3','.png'],'png');


vx=[];
vy=[];
figure(6)
for j=1:size(newC,1)
    xypoint=newC{j,1};
    for jj=1:size(xypoint,1)
        vx=[vx [newV(xypoint(jj,1),1); newV(xypoint(jj,2),1)]];
        vy=[vy [newV(xypoint(jj,1),2); newV(xypoint(jj,2),2)]];
    end
end

plot(x,y,'r+',vx,vy,'b-');

xlimit=[min(x)-20 max(x)+20];
ylimit=([min(y)-20 max(y)+20]);

xlim(xlimit);
ylim(ylimit);

saveas(gcf,[param.pathR,'final','.fig'],'fig');
saveas(gcf,[param.pathR,'final','.png'],'png');




