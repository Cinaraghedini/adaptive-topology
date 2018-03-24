% filename: eLocalN.m
% Purpose:  computes the local efficiency of a graph which is the average
% of each node local efficiency
% Input: 
% - mat - adjacency matrix
% - mwI - distance matrix
% - Options - graph options
% Output: 
% - eLocalG - array with the local efficiency of each node into the network
% Reference :
% Phys Rev Lett. 2001 Nov 5;87(19):198701. Epub 2001 Oct 17.
% Efficient behavior of small-world networks.
% Latora V1, Marchiori M.


function eLocalG = eLocalN(mat,mw,options)

if options.unweightedDeg
   [degG,~,~]=degrees(mat);
else
   [degG,~,~]=degrees_mod(mw);   
end
   
eLocalG = zeros(1,size(mw,1));

for v = 1:size(mw,1)  % for each node in the network
    
    nlist = kneighbors(mat,v,1); % list of node v neighbors
    sumE = 0;
    if (length(nlist) > 1) % if there is more than one neighboor
        for u=nlist
            for i=nlist
                 if mat(u,i)  %if there is an edge between neigbhors of v                    
                    sumE = 1/(mw(u,i)) + sumE;  
                end
            end
        end
    
        if options.unweightedDeg % weight of connected nodes
           coef = length(nlist)*(length(nlist)-1); 
        else
           coef = degG(v)*(length(nlist)-1); 
        end
        eLocalG(1,v) = (sumE/2)/coef;             
    end
    
end
