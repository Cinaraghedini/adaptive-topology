function properties = compute_propertiesGraph(g,options,param)


m=matrix(g);

distance = dist(g);

[m2 mw mwI] = initialize_matrix(distance,param,options);

matE=double(m);

matE=matE.*distance;  %matrix of weights 

properties = zeros(1,12);

if nv(g) > 1
    
    properties(1,1) = eGlobal(matE,distance,1);
    
    properties(1,2) = eLocal(m,matE,options);
    
    properties(1,3) = clusterCoefficientGraph(g);
    
    properties(1,4) = sum(weighted_clust_coeff(m))/nv(g);
    
    properties(1,5) = algebraic_connectivity_New(m,param.normalized); %not normalized = 0, normalized=1
    
    properties(1,6) = (max(subgrafosSize(g))/nv(g));
    
    properties(1,7) = sum(deg(g))/nv(g);
    
    properties(1,8:9) = robustnessMatrix(m,mw,param,'BC'); %network robustness   

    [data] = create_Ninf_Matrix(m);
    
    properties(1,10:11)= compute_numberPath_2(data,param);    

    properties(1,12) = robustnessMatrix_most(m,param); %network robustness

end