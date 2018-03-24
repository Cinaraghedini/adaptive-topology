% protocol to evaluate the performance of the combined control strategy
% that aims at improving the network robustness. The setModel file is the
% setup to the simulation environment. 

close all;
clear all;

clc;

% model parameter setup

graph_init;   % graph library
graphOptions; % graph parameters
setModel;    % model parameters

if isempty(param.networkList)
   param.networkList=transpose([1:param.numberNetworks].');
end

param.numberNetworks=length(param.networkList);

propertiesF =[]; % persists the average of the network evolution for each gain combination 

%--------------------------original code---------------------------

discT=[]; %persists the average of the disconnection iteration for each gain combination 

auxObstacles=obstacle_generator(param);

[param.failureT] = failureTime_generator(param); 

param.save=1;

[lg cg]=size(param.gain);

param.labels=[];
for i=1:size(param.gain,1)
    param.labels=[param.labels {['\sigma=' num2str(param.gain(i,1)) ',\psi=' num2str(param.gain(i,2))]}];
end

for k=3:3 %length(param.gain)

    discIt=[];
    
    properties=zeros(param.iteration,10);

    param.gainConnectivityController=param.gain(k,1);

    param.gainFormationController=param.gain(k,2);
    
    labelP=strcat(num2str(param.gainConnectivityController),'_',num2str(param.gainFormationController));
    
    param.pathR = [param.mainPath labelP '\'];

    param.fileName = labelP; 
    
    param.legend=[param.legend,param.labels(1,k)];
    
    param.XTickLabel=[param.XTickLabel,{labelP}];
    
     for ii=1:1%param.networkList
      
        param.obstacles=auxObstacles;
        
        param.listNodeFailure=[];
        
        param.network=ii;
        
        load([param.path num2str(param.network) '\position']);
    
        labels = num2str([1:size(position,1)]');
        x=transpose(position(:,1));
        y=transpose(position(:,2));
        
        figure(30);
        
        x=transpose(x);
        y=transpose(y);
        
        
        [vx,vy] = voronoi(x,y);
        
        plot(x,y,'r+',vx,vy,'b-');
        
        xlimit=[min(x)-3 max(x)+3];
        ylimit=([min(y)-3 max(y)+3]);
        
        xlim(xlimit);
        ylim(ylimit);

        
%        K=convhull(position(:,1), position(:,2));
        
        % Turn the hull into X,Y boundary values
        
%         BX = x(K); 
%         BY = y(K);
%         
%         % Do the special Voronoi limited by the data
%         [V, C, XY] = VoronoiLimit(x, y, 'bs_ext', [BX, BY], 'plot', 'on');
        
        
        distance = squareform(pdist(position,'euclidean'));
        
        [m mw mwI] = initialize_matrix(distance,param,options);
        
        figure(31);
        
        gplotVideo(m,position,labels);
        
        text(position(:,1),position(:,2), labels, 'VerticalAlignment','bottom','HorizontalAlignment','right','FontName', 'Segoe UI Semibold','Fontsize',15,'FontWeight','bold');
        
        newPosition=position;
        
         
        for i=1:size(position,1)

            listNi=kneighbors(m,i,1);
            positionNi=[];
            for j=listNi
                positionNi=[positionNi;position(j,:)];
            end 
            
           % [centroid, centroid_2] = centroid_computation_v7%(position(i,:),positionNi);
            
           % [centroid] = centroid_computation_v8(position(i,:),positionNi,param);
            
           [V, C, XY] = VoronoiLimit(x, y, 'bs_ext', [xmin ymin; xmin ymax; xmax ymax; xmax ymin], 'plot', 'on');
           
           
           
           
           
            if isempty(centroid)
                disp('No centroid found')
            else
                centroidPos(i,:)=centroid;
                
                newPosition(i,:)=centroidPos(i,:);
                
                xi=transpose(positionNi(:,1));
                yi=transpose(positionNi(:,2));
                
                figure(i);
                
                [vxi,vyi] = voronoi(xi,yi);
                
                plot(xi,yi,'r+',vxi,vyi,'b-');
                
                hold all;
                
                xlimit=[min(xi)-3 max(xi)+3];
                ylimit=([min(yi)-3 max(yi)+3]);
                
                xlim(xlimit);
                ylim(ylimit);
                              
                plot(centroidPos(i,1),centroidPos(i,2),'go');
                
                figure(i+param.networkSize)
                
                plot(x,y,'r+',vx,vy,'b-');
                
                xlimit=[min(x)-3 max(x)+3];
                ylimit=([min(y)-3 max(y)+3]);
                
                xlim(xlimit);
                ylim(ylimit);
                
                hold all;
                
                plot(centroidPos(i,1),centroidPos(i,2),'go');
                
            end
        end
     end
     
     figure(101);
     
     gplotVideo(m,newPosition,labels);
     
     text(newPosition(:,1),newPosition(:,2), labels, 'VerticalAlignment','bottom','HorizontalAlignment','right','FontName', 'Segoe UI Semibold','Fontsize',15,'FontWeight','bold');
     
     x=transpose(newPosition(:,1));
     y=transpose(newPosition(:,2));
     
     
     figure(200);
     
     [vx,vy] = voronoi(x,y);
     
     plot(x,y,'r+',vx,vy,'b-');
     
     xlimit=[min(x)-3 max(x)+3];
     ylimit=([min(y)-3 max(y)+3]);
     
     xlim(xlimit);
     ylim(ylimit);
     
end