function plotAnimation(step, area, time, qtyOfRobots, field)

    global subplotDim 
    
    plotScenario([area, time]);
    
    title(['Step ' num2str(step), ' - Area ' sprintf('%.3f', area), ' in ',sprintf('%.0f', 1000*time), ' miliseconds for ' sprintf('%0.f robots', qtyOfRobots)]);    
    
    if(nargin == 1 || ~exist('field', 'var') || isempty(field))
        field = getField();
    end
    
    if(~isempty(subplotDim))
        subplot(subplotDim(1), subplotDim(2), subplotDim(3));
    end
        
    axis(field);
    axis('square');
    color('b');
    grid on;
    hold on

    
        

end