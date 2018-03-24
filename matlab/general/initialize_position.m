function [position] = initialize_position(set) 

while true
    position = [((set.lowerB - set.upperB).*rand(set.networkSize,1) + set.upperB)  ((set.lowerB - set.upperB).*rand(set.networkSize,1) + set.upperB)];
    if (size((unique(strcat(num2str(position(:,1)),num2str(position(:,2))),'rows')),1) == set.networkSize)
        return;
    end
end