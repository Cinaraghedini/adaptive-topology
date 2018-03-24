function XY = createPointsUniform(position,param)


xlimUpper=max(position(:,1));
xlimLower=min(position(:,1));

ylimUpper=max(position(:,2));
ylimLower=min(position(:,2));

if xlimLower<0 && xlimUpper<0
    segmentX=abs(abs(xlimUpper)-(abs(xlimLower)));
else
    if xlimLower<0 && xlimUpper>0
        segmentX=abs(abs(xlimUpper)+(abs(xlimLower)));
    else
        segmentX =xlimUpper-xlimLower;
    end
end
        
if ylimLower<0 && ylimUpper<0
    segmentY=abs(abs(ylimUpper)-(abs(ylimLower)));
else
    if ylimLower<0 && ylimUpper>0
        segmentY=abs(abs(ylimUpper)+(abs(ylimLower)));
    else
        segmentY =ylimUpper-ylimLower;
    end
end

segmentX=round(segmentX);
segmentY=round(segmentY);

area=segmentX*segmentY;
points= round(area/size(position,1));

v=zeros(1,round(area));

pointV = [1 + (0:size(position,1)-1)*points]';

v(pointV)=1;

m=reshape(v,segmentX,segmentY);

[l c]=find(m);

XY=[l c];