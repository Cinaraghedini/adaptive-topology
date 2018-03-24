function [V C]=voronoiBounded(xy,param)

[V, C] = voronoin(xy); 

Vbounded=V;
Cbounded=C;

outside_Points=[];

for i=1:size(C,1)
    pos_i=xy(i,:);
    P=C{i,1};
    for ii=P
        posii=V(ii,:);
        if isinf(posii)
        end 
        distancei_ii=pdist([posii; pos_i],'euclidean');
        if distancei_ii > param.range
           new_ij(1,1) = pos_i(1,1) + ((posii(1,1) - pos_i(1,1))*param.range)/distancei_ii;
           new_ij(1,2) = pos_i(1,2) + ((posii(1,2) - pos_i(1,2))*param.range)/distancei_ii;
           new= pos_i+((posii-pos_i)*param.range)/distancei_ii;
        end
    end
end