% filename: plotCoverage_Integrated.m
% Purpose:  plot and persist the coverage area evolution for all gain
% setting combinations
% Input: 
% - data - matrix containing the coverage data evoluation
% - param - parametrization struct

function plotCoverage_Integrated(coverageData,param)

close all;

set(0,'defaultfigurecolor',[1 1 1])

strTitle=sprintf('N=%d, Area=$%d^2$, Iterations=%d', param.networkSize,param.area,param.numberNetworks);


figure(1);   %density

for i=1:size(coverageData,1)
    
    data=coverageData{i,1};
 
    plot(data(:,1),data(:,2),'LineStyle','-','LineWidth',2,'Marker','x');
    
    hold all;
end

title(strTitle,'Interpreter','latex','FontSize',14);

h = legend(param.legend,'Location','Best','FontSize',14);

ylabel('Density','FontSize',14);
xlabel(param.labelX,'FontSize',14);
ylim([0 1]); 

gName=[param.mainPath,'graphic_Density_',param.fileName,'_',num2str(param.networkSize)];

saveas(gcf,[gName,'.fig'],'fig');
export_fig(gName,'-eps','-png');


figure(2);  % Area coverage rate 
 
for i=1:size(coverageData,1)
    
    data=coverageData{i,1};
    plot(data(:,1),data(:,3),'LineStyle','-','LineWidth',2,'Marker','x');
    hold all;
end

title(strTitle,'Interpreter','latex','FontSize',14);

h = legend(param.legend,'Location','Best','FontSize',14);

ylabel('Area coverage rate','FontSize',14);

xlabel(param.labelX,'FontSize',14);

gName=[param.mainPath,'graphic_areaCovered_',param.fileName,'_',num2str(param.networkSize)];

saveas(gcf,[gName,'.fig'],'fig');
export_fig(gName,'-eps','-png');


figure(3);   % Averaged area coverage rate 

maxR=[];
for i=1:size(coverageData,1)
    data=coverageData{i,1};
    maxR=[maxR; data(:,4)];
    plot(data(:,1),data(:,4),'LineStyle','-','LineWidth',2,'Marker','x');
    hold all;

end

title(strTitle,'Interpreter','latex','FontSize',14);

h = legend(param.legend,'Location','Best','FontSize',14);

ylabel('Area coverage rate - averaged','FontSize',14);

xlabel(param.labelX,'FontSize',14);
ylim([0 max(maxR)+1]);

gName=[param.mainPath,'graphic_averagedArea_',param.fileName,'_',num2str(param.networkSize)];

saveas(gcf,[gName,'.fig'],'fig');
export_fig(gName,'-eps','-png');


figure(4);   % average distance to centroid
 
maxR=[];
for i=1:size(coverageData,1)
    data=coverageData{i,1};
    maxR=[maxR; data(:,5)];
    plot(data(:,1),data(:,5),'LineStyle','-','LineWidth',2,'Marker','x');
    hold all;
end

title(strTitle,'Interpreter','latex','FontSize',14);

h = legend(param.legend,'Location','Best','FontSize',14);

ylabel('$avgDist(\mathcal{G})$','Interpreter','latex','FontSize',14);

xlabel(param.labelX,'FontSize',14);

ylim([0 max(maxR)+1]);

gName=[param.mainPath,'graphic_AVGCentroid_',param.fileName,'_',num2str(param.networkSize)];

saveas(gcf,[gName,'.fig'],'fig');
export_fig(gName,'-eps','-png');


figure(5);   % fraction of covered points 

for i=1:size(coverageData,1)
    data=coverageData{i,1};
    plot(data(:,1),data(:,6),'LineStyle','-','LineWidth',2,'Marker','x');
    hold all;
end

title(strTitle,'Interpreter','latex','FontSize',14);

h = legend(param.legend,'Location','Best','FontSize',14);

ylabel('$f_{Spoints_c}$','Interpreter','latex','FontSize',16);
xlabel(param.labelX,'FontSize',14);
ylim([0 1]); 



gName=[param.mainPath,'graphic_Fraction_Inside_',param.fileName,'_',num2str(param.networkSize)];

saveas(gcf,[gName,'.fig'],'fig');
export_fig(gName,'-eps','-png');


figure(6);    % fraction of uncovered points 

for i=1:size(coverageData,1)
    data=coverageData{i,1};
    plot(data(:,1),data(:,8),'LineStyle','-','LineWidth',2,'Marker','x');
    hold all;
end
 
title(strTitle,'Interpreter','latex','FontSize',14);

h = legend(param.legend,'Location','Best','FontSize',14);

ylabel('$f_{Spoints_u}$','Interpreter','latex','FontSize',16);

xlabel(param.labelX,'FontSize',14);
ylim([0 1]); 

gName=[param.mainPath,'graphic_Fraction_Outside_',param.fileName,'_',num2str(param.networkSize)];

saveas(gcf,[gName,'.fig'],'fig');
export_fig(gName,'-eps','-png');


figure(7);   % fraction of overlap points 

for i=1:size(coverageData,1)
    data=coverageData{i,1};
    plot(data(:,1),data(:,10),'LineStyle','-','LineWidth',2,'Marker','x');
    hold all;
end

title(strTitle,'Interpreter','latex','FontSize',14);

h = legend(param.legend,'Location','Best','FontSize',14);

ylabel('$Spoints_{overlap}$','Interpreter','latex','FontSize',16);

xlabel(param.labelX,'FontSize',14);
ylim([0 1]); 

gName=[param.mainPath,'graphic_Fraction_Overlap_',param.fileName,'_',num2str(param.networkSize)];

saveas(gcf,[gName,'.fig'],'fig');
export_fig(gName,'-eps','-png');

close all;