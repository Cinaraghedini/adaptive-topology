function ret = getUnique(vector)

    vector(isnan(vector)) = [];
    ret = unique(vector);

end