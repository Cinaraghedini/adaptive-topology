function point = getPoint(x, y, robotOrig, robotDest, polygonOrder)    
    point = nan(5, 1);
    
    point(1) = x;
    point(2) = y;
    point(3) = robotOrig;
    point(4) = robotDest;    
    point(5) = polygonOrder;
    
    point = point';
    
end