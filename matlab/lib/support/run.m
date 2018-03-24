function run
       
    clear -global subplotDim exit freeze
    close all
    averageSpeed = 10; %m/s
    field = [-150, 200, -150, 200];
    
    exit = false;
    freeze = false;
    
    qtyOfScenes = 100;
    
    R = 20;
    scenario = 1;
    load(strcat('database/',num2str(scenario),'/position.mat'))
    robots = [position repmat(R, length(position), 1)];
    
    animation(robots, averageSpeed, field, qtyOfScenes);
end


function animation(robotsXYR, averageSpeed, field, qtyOfScenes)

    clear -global exit freeze scenes

    [m n] = size(robotsXYR);
    
    h = figure('units','normalized','outerposition',[0 0 1 1]);
    set(h,'KeyPressFcn',@callbackFcn);
    
    global subplotDim exit freeze scenes
    
    
    
    scenes = java.util.HashMap;
    scenes.put(0, robotsXYR);
    
    for(i = 1 : qtyOfScenes)
        
        if(exit)
            break;
        end
        
        clf
%         disp(['>> Passo ' num2str(i) ':'])
        [area time] = calcArea(robotsXYR);
        angles = 2 * pi * rand(m,1);
        
        robotsXYR(:,1) = robotsXYR(:,1) + averageSpeed * cos(angles);
        robotsXYR(:,2) = robotsXYR(:,2) + averageSpeed * sin(angles);
        
        scenes.put(i, robotsXYR);
               
        plotAnimation(i, [area time], field);
        drawnow
        
        while(freeze)
            pause(0.01);
        end       
        
    end
    
end

function callbackFcn(h,evt)
    global exit freeze scenes
    
    if(strcmp(evt.Key, 's'))
        disp('Dados salvos em animation.mat')
        save('./animation.mat', 'scenes');

        return;
    end
    
    if(strcmp(evt.Key, 'escape'))
        exit = true;
        close all;
        clc;        
        return;
    end   
        
    if(freeze)
        freeze = false;
        else
            freeze = true;
    end
end