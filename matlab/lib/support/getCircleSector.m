function [X, Y] = getCircleSector(center, radius, fromAngle, toAngle)

    x0 = center(1);
    y0 = center(2);
        
    delta = toAngle - fromAngle;
%     toAngle = fromAngle + delta;
    
    X = x0;
    Y = y0;    
    
    for i = fromAngle : delta / 100 : toAngle
        X = [X, x0 + radius * cos(i * pi /180)]; %#ok<*AGROW>
        Y = [Y, y0 + radius * sin(i * pi /180)];
    end
    
     X = [X, x0 + radius * cos(toAngle * pi /180)];
     Y = [Y, y0 + radius * sin(toAngle * pi /180)];
        
end