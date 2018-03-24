function plotScenario(outcome, field)

    global inner outer outerPolygons innerPolygons sectors robots overlapedInnerPolygons keep
    
    if(isempty(keep))
        keep = false;
    end
    
    if(~keep)
        figure('units','normalized','outerposition',[0 0 0.5 1])
    end

    if(nargin == 2)
        plotGrid(robots, field);
    else
        plotGrid(robots);
    end
    
    [qtyOfRobots, ~] = size(robots);
         
    title(['Area ' sprintf('%.2f', outcome(1)), ' in ',sprintf('%.2f', 1000*outcome(2)), ' miliseconds for ' sprintf('%d robots', qtyOfRobots)]);    

    if(~isempty(outer))
        hold on
        plot(outer(:, 1), outer(:, 2), 'r*');
        plot(inner(:, 1), inner(:, 2), 'b.');

        plotPolygons(outerPolygons, 'r');
        plotPolygons(innerPolygons, 'b');
        plotPolygons(overlapedInnerPolygons, 'k', true);

        plotSectors(sectors);
    end
    
    hold off
    grid minor

end