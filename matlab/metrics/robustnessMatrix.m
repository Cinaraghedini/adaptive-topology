function [data]=robustnessMatrix(matrix,matrixW,param,type)

[nd nl]=size(matrix);

if param.fractionIteration
   nrNodes=1;
else
   nrNodes=round(param.fractionIteration*nd);
end

data=zeros(1,2);

if nrNodes<nd
    
    count=0;
    
    g=graph;
    
    connected=true;
    
    while connected && nd>nrNodes
        
        set_matrix(g,matrix);
        connected=isconnected(g);
        nd=nv(g);
        if connected
            
            count=count+1;

            if strcmp(type,'Failures')
                idx=randi([1 nd]);
            else
                if strcmp(type,'BC')
                    BC=[betweenness_wei(matrixW)  [1:nd].'];   % compute the node's BC
                    [ranking pos] = sort(BC(:,1),'descend');  % ordering nodes according to their BC ranking
                    idx=pos(1:nrNodes);                    
                end
            end
            
            matrix(:,idx)=[];
            matrix(idx,:)=[];
            matrixW(:,idx)=[];
            matrixW(idx,:)=[];
            
        end
 
    end
    
    if nv(g)==0
        fnv=1;
    else
        fnv=max(subgrafosSize(g))/nv(g);
    end
    data=[count/nl fnv];
    free(g);
    
end
