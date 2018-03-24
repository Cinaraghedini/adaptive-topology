function singlePath_nodes = find_2hneighbor(data,param)

singlePath_nodes=[];

for i=1:size(data,1)
    
    fractionPath=0;
    
    node=data{i,1}{1,1};
     
    nNeighbors=size(data{i,1}{1,2},2);
    
    dds=[];
     
    if nNeighbors > 0
    
        for j=2:nNeighbors+1  % neighbors
           if ~isempty(data{i,1}(j,2))
              dds=[dds; cell2mat(data{i,1}(j,2))];
           end
        end
            
    end
    
    countNeigh=0;
 
    if ~isempty(dds)
        
        idNodes=transpose(unique(dds(:,1)));
         
        for j=idNodes
            if isempty(find(data{i,1}{1,2}==j)) && j~=node % se não é vizinho direto do nodo
                nodesPath=find(dds(:,1)==j);
                if size(nodesPath,1) == param.nrPaths
                    singlePath_nodes=[singlePath_nodes;i  j];
                end
            end
        end
        
    end
end
        