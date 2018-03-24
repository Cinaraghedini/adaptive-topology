function [data]=robustnessMatrix_most(matrix,param)

[nd nl]=size(matrix);

if param.fractionIteration
   nrNodes=1;
else
   nrNodes=round(param.fractionIteration*nd);
end

dataA=[];

if nrNodes<nd
    
    for i=1:size(matrix,1)
        
        v=i;

        g=graph;
    
        set_matrix(g,matrix);

        delete(g,v);
                
        dataA=[dataA; [v max(subgrafosSize(g))/nv(g)]];
        
        free(g);
        
    end
      
end

[S I]=sort(dataA(:,2));

data=S(1,1);


