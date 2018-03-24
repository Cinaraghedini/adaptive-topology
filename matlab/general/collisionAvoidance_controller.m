function [output] = collision_avoidance(position,param)

positionAvoidance=[position; param.obstacles]; 

distance = squareform(pdist(positionAvoidance,'euclidean'));

output=zeros(size(position,1),2);  %nv número de robôs

for i=1:size(position,1)
    pos=find(distance(i,:)<=param.minDistObstacle);
    for j=pos
        K=0;
        if j~=i
            K = distance(i,j)^2 - ((param.minDistObstacle^3)/distance(i,j));
        end
        output(i,1) = output(i,1) - K * (positionAvoidance(i,1)-positionAvoidance(j,1));
        output(i,2) = output(i,2) - K * (positionAvoidance(i,2)-positionAvoidance(j,2));
    end
end