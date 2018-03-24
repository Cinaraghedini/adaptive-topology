function t = lesserOrEqualThan(a, b)

    precision = 5;
    
    a = round(a*power(10, precision));
    b = round(b*power(10, precision));
    
    t = (a <= b);
    
end