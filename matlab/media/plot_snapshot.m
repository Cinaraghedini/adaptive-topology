% filename: plot_Video.m
% Purpose:  generating a figure from one timestamp data (network state)
% Input:
% - data: matrix containing the simulation time, node positions, state of node (vulnerable or not)
% - param: parametrization setup
% - theta: network robustness level for title
% - nodeToFail: contains the node that will fail, if any, at the next simulation iteration. It is used for highliting it before failure
% - nodeToObstacle: contains the node that failed in the previous
% iterations, if any. It is used for highliting it and showing it as an obstacle.
% Output:
% - currentFrame: figure showing the network status with all the information about the node state


function currentFrame = plot_snapshot(data,theta,nodeToFail,node2obstacle,param)

warning('off','all')

position=data(:,3:4); %node positions

[matrix] = initialize_matrixAdj(position,param); % adjacency matrix initialization

if param.frameOn
    set(0,'DefaultFigureVisible','on');
else
    set(0,'DefaultFigureVisible','off');
end

legendDsc=[];

% title composition

headerString=strcat('(\sigma= ',num2str(param.gain(param.ref,1)),',\psi= ',num2str(param.gain(param.ref,2)), ',\zeta=', num2str(param.gain(param.ref,3)),')    \Theta= ',sprintf('%1.3f',theta),'    t= ',sprintf('%1.5f',(data(1,1)+param.idx)));

plot(param.obstacles(:,1),param.obstacles(:,2),'LineStyle','none','MarkerFaceColor','w','MarkerEdgeColor','w','Marker','.','MarkerSize',0.1);

% plotting obstacles

if param.showObstacles
    text(param.obstacles(:,1),param.obstacles(:,2),'I', 'FontName','DFKai-SB','Fontsize',7,'FontWeight','bold');
    legendDsc=[legendDsc {'| obst'}];    
end

hold all;

% plotting the graph specified by matrix and x y coordinates at position.
gplotVideo(matrix,position,'');  

legendDsc=[legendDsc {'links'}];

hold all;

% highlighting node that will fail in the next iteration, if any.
if ~isempty(node2obstacle)
    plot(node2obstacle(:,1),node2obstacle(:,2),'LineStyle','none','MarkerFaceColor',[0.800000011920929 0.600000023841858 0],'MarkerEdgeColor',[0.800000011920929 0.600000023841858 0],'MarkerSize',3,'Marker','o');
    legendDsc=[legendDsc {'failed node'}];
    hold all;
end


if param.color % if vulnerable nodes should be highlighted
    
    [ddsplot] = findPos_State(data,1); % finds vulnerable nodes
    
    if ~isempty(ddsplot) % if there are vulnerable nodes
        plot(ddsplot(:,3),ddsplot(:,4),'LineStyle','none','MarkerFaceColor',[1 0.200000002980232 0.200000002980232],'MarkerEdgeColor',[1 0.200000002980232 0.200000002980232],'MarkerSize',3,'Marker','o');
        legendDsc=[legendDsc {'vulnerable nodes'}];
        hold all;
        
    end
end
% else
[ddsplot] = findPos_State(data,0); % finds non vulnerable nodes
plot(ddsplot(:,3),ddsplot(:,4),'LineStyle','none','MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 1],'MarkerSize',3,'Marker','o');
hold all;

legendDsc=[legendDsc {'nodes'}];


title(headerString, 'FontSize',14);

h=legend(legendDsc{1,:},'Orientation','horizontal','Position','Location','SouthOutside');

xlim([0 500]);
ylim([-10 60]);

if ~isempty(nodeToFail)
    plot(nodeToFail(1,1),nodeToFail(1,2),'MarkerEdgeColor',[0 1 0],'MarkerSize',10,...
        'Marker','o',...
        'LineWidth',2,...
        'LineStyle',':');
end

set(h,'FontSize',10);

legend boxoff;

if param.frameOn
    currentFrame = getframe(gcf);
else
    currentFrame = hardcopy(gcf, '-dzbuffer', '-r0');
end