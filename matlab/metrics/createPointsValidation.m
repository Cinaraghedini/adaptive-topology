% filename: createPointsValidation.m
% Purpose:  creates uniformly distributed random spots considering a square area around the nodes positions 
% the number of spots to be generated is given by param.nrPoints
% Input: 
% - position - x and y coordinates for each network agent
% - param - parametrization struct
% Output: 
% XY - (n,2) matrix with the x,y position for each of the n spots 


function XY = createPointsValidation(position,param)

xlimUpper=max(position(:,1));
xlimLower=min(position(:,1));

ylimUpper=max(position(:,2));
ylimLower=min(position(:,2));

xLower=xlimLower-param.range;
yLower=ylimLower-+param.range;
xUpper=xlimUpper+param.range;
yUpper=ylimUpper+param.range;

if xLower<0 && xUpper<0
    segmentX=abs(abs(xUpper)-(abs(xLower)));
else
    if xLower<0 && xUpper>0
        segmentX=abs(abs(xUpper)+(abs(xLower)));
    else
        segmentX =xUpper-xLower;
    end
end
        
if yLower<0 && yUpper<0
    segmentY=abs(abs(yUpper)-(abs(yLower)));
else
    if yLower<0 && yUpper>0
        segmentY=abs(abs(yUpper)+(abs(yLower)));
    else
        segmentY =yUpper-yLower;
    end
end
        
xPoints= segmentX/param.nrPoints;
yPoints= segmentY/param.nrPoints;

X = [xLower + (0:param.nrPoints)*xPoints]';
Y = [yLower + (0:param.nrPoints)*yPoints]';

XY=[];
n=size(Y);
for i=1:size(X,1)
    r=repmat(X(i,1),1,n)';
    XY=[XY; [r Y]];
end