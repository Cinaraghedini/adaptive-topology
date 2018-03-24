
function [ param ] = voronoiCentroid(position,param,options)


labels = num2str([1:size(position,1)]');

x=position(:,1);
y=position(:,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% original voronoi cells

figure(1);

[vx,vy] = voronoi(x,y);

plot(x,y,'b*');

xlimit=[min(x)-10 max(x)+10];
ylimit=([min(y)-10 max(y)+10]);

xlim(xlimit);
ylim(ylimit);
    
[centroids]=define_centroid_local(position,options,param);

hold all;

%plot(centroids(:,1),centroids(:,2),'ro');

%hold all;

plot(vx,vy,'k-');

%legend('robot','centroid','voronoi Tesselation')

legend('robot','voronoi Tesselation')

%saveas(gcf,[param.pathR,'originalCentroid_',num2str(param.network),'.png'],'png');
%saveas(gcf,[param.pathR,'originalCentroid_',num2str(param.network),'.eps'],'psc2');
export_fig([param.pathR,'centroid_',num2str(param.network)], '-eps','-png');

figure(2);

[centroids]=define_centroid_local(position,options,param);

plot(centroids(:,1),centroids(:,2),'b*');

hold all;

plot(vx,vy,'k-');

xlim(xlimit);
ylim(ylimit);

legend('robot','voronoi Tesselation')

%saveas(gcf,[param.pathR,'centroid_',num2str(param.network),'.png'],'png');

%saveas(gcf,[param.pathR,'centroid_',num2str(param.network),'.eps'],'psc2');

export_fig([param.pathR,'originalCentroid_',num2str(param.network)], '-eps','-png');


for i=1:4

    x=centroids(:,1);
    y=centroids(:,2);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % centroidal voronoi cells

    [centroids]=define_centroid_local(centroids,options,param);
    
    figure(i+3)

    [vx,vy] = voronoi(x,y);

    plot(x,y,'b*',vx,vy,'k-');

    xlim(xlimit);
    ylim(ylimit);

    %saveas(gcf,[param.pathR,'centroidalDiagram_',num2str(param.network),'_',num2str(i),'.png'],'png');
    
    %saveas(gcf,[param.pathR,'centroidalDiagram_',num2str(param.network),'_',num2str(i),'.eps'],'psc2');

    export_fig([param.pathR,'centroidalDiagram_',num2str(param.network),'_',num2str(i)], '-eps','-png');


    
end
