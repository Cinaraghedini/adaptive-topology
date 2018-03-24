function equal  = isEqualTo(valueA, valueB)
    
    precision = 5;
    
    valueA = round(valueA * power(10, precision));
    valueB = round(valueB * power(10, precision));
    
    equal = valueA == valueB;

end