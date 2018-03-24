function evalComplexity(quantitiesOfRobots, rangeOfEachRobot)
    
    close all    
    clear -global
    
    m = length(quantitiesOfRobots);
    n = max(quantitiesOfRobots);
    
    timesInSecs = nan(m, 1);
    areas = nan(m, 1);
    
    for(i = 1 : m)
        clear robots        
        fieldLength = 1.5 * quantitiesOfRobots(i);
        
        robots = fieldLength * rand(n, 2);
        robots = [robots repmat(rangeOfEachRobot, n, 1)];
        
        samples = robots(1 : quantitiesOfRobots(i), :);
        tic;
        area = calcArea(samples);
        timeInSecs = toc;
        timesInSecs(i) = timeInSecs;
        areas(i) = area;
        
        disp(sprintf('Área de %d robôs: %d em %d segundos.',  quantitiesOfRobots(i), area, timesInSecs(i)));
        
    end
         
    plot(quantitiesOfRobots, timesInSecs, 'bo-');
    title('Tempo de Execução em s');
    xlabel('Quantidade de robôs');
    ylabel('Tempo (s)');
    hold on
    
    px = linspace(min(quantitiesOfRobots), max(quantitiesOfRobots), 10000);
    
%     p1 = polyfit(quantitiesOfRobots',timesInSecs, 1);
    p2 = polyfit(quantitiesOfRobots',timesInSecs, 2);
    p3 = polyfit(quantitiesOfRobots',timesInSecs, 3);
    
    py = polyval(p2, px);
    py2 = polyval(p2, px) + 10;
    py3 = polyval(p3, px) + 20;
         
    plot(px, py2, 'r--');    
    
    grid on   
    save('./saves/complexity.mat', 'quantitiesOfRobots', 'timesInSecs', 'areas', 'robots', 'p2');

end