function field = getField(robots)

    r = max(robots(:,3));
    
    x = [min(robots(:,1)) - r,  max(robots(:,1)) + r];
    y = [min(robots(:,2)) - r,  max(robots(:,2)) + r];
    
    range = abs(max(diff(x), diff(y)));

    field = [x(1), x(1)+range, y(1), y(1)+range];    
    
end