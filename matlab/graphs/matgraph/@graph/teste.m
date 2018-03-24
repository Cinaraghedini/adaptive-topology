
% ATENTAR DELECAO ARESTA Cc

t = cputime
graph_init
g = graph
n=5

iteracao=0


% criação grafo

add(g,1,2)
add(g,2,3)
add(g,3,4)
add(g,4,5)
add(g,2,5)
ndraw(g)

n = nv(g);   % número de vértices de g

q_init(n+1) 
track = zeros(1,n)
track(v) = v

path= zeros(n,n)
path(v,:) = v
path(v,:) = v


q_push(v) 
 
while(q_size > 0)
     t = q_pop_front
     if t==u
        break
     end
     push_list = neighbors(g,t)
     for s = push_list
        if (track(s) == 0)
            track(s) = t
            q_push(s)
        end
        if (track(s) ~= 0)
            path(t,s)=t
        end
    end
end
disp(track)
disp(path)
waitforbuttonpress
  

% Shortest Path Lenght
    
    %MedShortestPathInd=zeros(1,n)
    %L=0
    %for v=1:n
    %    SumShortestPathInd=0
    %    for u=1:n
    %        if (v ~= u)
    %            shortestPath = find_path(g,v,u)
    %            [linha coluna] = size(shortestPath)
    %            SumShortestPathInd =  SumShortestPathInd+(coluna-1)
    %            
    %        end
    %    end
    %    MedShortestPathInd(1,v)=SumShortestPathInd/(n-1) % Average Shortest Path Lenght -  Individual 
    %end
    %L=median(MedShortestPathInd) % Average Shortest Path Lenght -  Graph
