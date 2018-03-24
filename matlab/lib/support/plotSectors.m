function plotSectors(sectors)

    if(isempty(sectors))
        return;
    end

    sectorsIds = unique(sectors(:,3));
    
    for i = 1 : length(sectorsIds)
        sector = sectors(sectors(:, 3) == sectorsIds(i), 1:2);
        x0 = sector(1, 1);
        y0 = sector(1, 2);
        sector(1,:) = [];
        patch([x0 sector(:,1)' x0], [y0 sector(:,2)' y0], 'y', 'FaceAlpha',0.25,'EdgeColor','None');
    end
end