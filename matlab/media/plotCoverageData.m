% filename: plotCoverageData.m
% Purpose:  plot and persist the coverage area evolution for a specific
% gain configuration
% Input: 
% - data - matrix containing the coverage data evoluation
% - param - parametrization struct


function plotCoverageData(data,param)

close all;

% title 

strTitle=sprintf('N=%d, Area=$%d^2$, Iterations=%d', param.networkSize,param.area,param.numberNetworks);

figure(1);   % network density
 
plot(data(:,1),data(:,2),'LineStyle','-','LineWidth',2,'Marker','x');

title(strTitle,'Interpreter','latex','FontSize',14);

ylabel('Density','FontSize',14);
xlabel(param.labelX,'FontSize',14);
ylim([0 1]); 

annotation('textbox',[0.684928571428568 0.816666666666672 0.191071428571429 0.0690476190476191],...
    'String',strcat('Gains: \sigma=',sprintf('%d',param.gainConnectivityController),',\psi=',sprintf('%d',param.gainRobustnessControl),',\zeta=',sprintf('%d', param.gainCoverageController) ),...
    'FontSize',12,...
    'FontName','Cambria',...
    'FitBoxToText','off',...
    'LineStyle','none');

gName=[param.pathR,'graphic_Density_',param.fileName,'_',num2str(param.networkSize)];

saveas(gcf,[gName,'.fig'],'fig');
saveas(gcf,[gName,'.png'],'png');
saveas(gcf,[gName,'.eps'],'psc2');


figure(2);   % Area coverage rate
 
plot(data(:,1),data(:,3),'LineStyle','-','LineWidth',2,'Marker','x');

title(strTitle,'Interpreter','latex','FontSize',14);

ylabel('Area coverage rate','FontSize',14);
xlabel(param.labelX,'FontSize',14);

annotation('textbox',[0.684928571428568 0.816666666666672 0.191071428571429 0.0690476190476191],...
    'String',strcat('Gains: \sigma=',sprintf('%d',param.gainConnectivityController),',\psi=',sprintf('%d',param.gainRobustnessControl),',\zeta=',sprintf('%d', param.gainCoverageController) ),...
    'FontSize',12,...
    'FontName','Cambria',...
    'FitBoxToText','off',...
    'LineStyle','none');
gName=[param.pathR,'graphic_coveredArea_',param.fileName,'_',num2str(param.networkSize)];

saveas(gcf,[gName,'.fig'],'fig');
saveas(gcf,[gName,'.png'],'png');
saveas(gcf,[gName,'.eps'],'psc2');


figure(3);   % Averaged area coverage rate 
 
plot(data(:,1),data(:,4),'LineStyle','-','LineWidth',2,'Marker','x');

title(strTitle,'Interpreter','latex','FontSize',14);

ylim([0 max(data(:,4))+1]); 

ylabel('Area coverage rate - averaged','FontSize',14);
xlabel(param.labelX,'FontSize',14);

annotation('textbox',[0.684928571428568 0.816666666666672 0.191071428571429 0.0690476190476191],...
    'String',strcat('Gains: \sigma=',sprintf('%d',param.gainConnectivityController),',\psi=',sprintf('%d',param.gainRobustnessControl),',\zeta=',sprintf('%d', param.gainCoverageController) ),...
    'FontSize',12,...
    'FontName','Cambria',...
    'FitBoxToText','off',...
    'LineStyle','none');

gName=[param.pathR,'graphic_averageCoveredArea_',param.fileName,'_',num2str(param.networkSize)];

saveas(gcf,[gName,'.fig'],'fig');
saveas(gcf,[gName,'.png'],'png');
saveas(gcf,[gName,'.eps'],'psc2');

figure(4);   % average distance to centroid
 
plot(data(:,1),data(:,5),'LineStyle','-','LineWidth',2,'Marker','x');

title(strTitle,'Interpreter','latex','FontSize',14);

ylabel('$avgDist(\mathcal{G})$','Interpreter','latex','FontSize',14);

xlabel(param.labelX,'FontSize',14);

ylim([0 max(data(:,5)+1)]);

annotation('textbox',[0.684928571428568 0.816666666666672 0.191071428571429 0.0690476190476191],...
    'String',strcat('Gains: \sigma=',sprintf('%d',param.gainConnectivityController),',\psi=',sprintf('%d',param.gainRobustnessControl),',\zeta=',sprintf('%d', param.gainCoverageController) ),...
    'FontSize',12,...
    'FontName','Cambria',...
    'FitBoxToText','off',...
    'LineStyle','none');

gName=[param.pathR,'graphic_AVGCentroid_',param.fileName,'_',num2str(param.networkSize)];

saveas(gcf,[gName,'.fig'],'fig');
saveas(gcf,[gName,'.png'],'png');
saveas(gcf,[gName,'.eps'],'psc2');


figure(5);  % fraction of covered, uncovered and overlap points 
 
plot(data(:,1),data(:,6),'LineStyle','-','LineWidth',2,'Marker','x');
hold all;
plot(data(:,1),data(:,8),'LineStyle','-','LineWidth',2,'Marker','o');
hold all;
plot(data(:,1),data(:,10),'LineStyle','-','LineWidth',2,'Marker','*');

title(strTitle,'Interpreter','latex','FontSize',14);

h = legend('f_{Spoints_c}','f_{Spoints_u}','Spoints_{overlap}','Interpreter','latex','Location','Best');

xlabel(param.labelX,'FontSize',14);
ylabel('Fraction of points','FontSize',14);

ylim([0 1]); 

annotation('textbox',[0.684928571428568 0.816666666666672 0.191071428571429 0.0690476190476191],...
    'String',strcat('Gains: \sigma=',sprintf('%d',param.gainConnectivityController),',\psi=',sprintf('%d',param.gainRobustnessControl),',\zeta=',sprintf('%d', param.gainCoverageController) ),...
    'FontSize',12,...
    'FontName','Cambria',...
    'FitBoxToText','off',...
    'LineStyle','none');

gName=[param.pathR,'graphic_Fraction_',param.fileName,'_',num2str(param.networkSize)];

saveas(gcf,[gName,'.fig'],'fig');
saveas(gcf,[gName,'.png'],'png');
saveas(gcf,[gName,'.eps'],'psc2');


figure(6);   % fraction of normalized covered, uncovered and overlap points 
 
plot(data(:,1),data(:,7),'LineStyle','-','LineWidth',2,'Marker','x');
hold all;
plot(data(:,1),data(:,9),'LineStyle','-','LineWidth',2,'Marker','o');
hold all;
plot(data(:,1),data(:,10),'LineStyle','-','LineWidth',2,'Marker','*');

title(strTitle,'Interpreter','latex','FontSize',14);

h = legend('f_{Spoints_c}','f_{Spoints_u}','Spoints_{overlap}','Interpreter','latex','Location','Best');

ylabel('Fractions of points - Normalized','FontSize',14);
xlabel(param.labelX,'FontSize',14);
ylim([0 1]); 

annotation('textbox',[0.684928571428568 0.816666666666672 0.191071428571429 0.0690476190476191],...
    'String',strcat('Gains: \sigma=',sprintf('%d',param.gainConnectivityController),',\psi=',sprintf('%d',param.gainRobustnessControl),',\zeta=',sprintf('%d', param.gainCoverageController) ),...
    'FontSize',12,...
    'FontName','Cambria',...
    'FitBoxToText','off',...
    'LineStyle','none');


gName=[param.pathR,'graphic_Fraction_Norm_',param.fileName,'_',num2str(param.networkSize)];

saveas(gcf,[gName,'.fig'],'fig');
saveas(gcf,[gName,'.png'],'png');
saveas(gcf,[gName,'.eps'],'psc2');



figure(7);  % fraction of covered, uncovered and overlap points regarding the size of the giant component
 
plot(data(:,1),data(:,11),'LineStyle','-','LineWidth',2,'Marker','x');
hold all;
plot(data(:,1),data(:,13),'LineStyle','-','LineWidth',2,'Marker','o');
hold all;
plot(data(:,1),data(:,15),'LineStyle','-','LineWidth',2,'Marker','*');

title(strTitle,'Interpreter','latex','FontSize',14);

h = legend('f_{Spoints_c}','f_{Spoints_u}','Spoints_{overlap}','Interpreter','latex','Location','Best');

ylabel('Fraction S(g)','FontSize',14);
xlabel(param.labelX,'FontSize',14);
ylim([0 1]); 

annotation('textbox',[0.684928571428568 0.816666666666672 0.191071428571429 0.0690476190476191],...
    'String',strcat('Gains: \sigma=',sprintf('%d',param.gainConnectivityController),',\psi=',sprintf('%d',param.gainRobustnessControl),',\zeta=',sprintf('%d', param.gainCoverageController) ),...
    'FontSize',12,...
    'FontName','Cambria',...
    'FitBoxToText','off',...
    'LineStyle','none');

gName=[param.pathR,'graphic_Fraction_S_',param.fileName,'_',num2str(param.networkSize)];

saveas(gcf,[gName,'.fig'],'fig');
saveas(gcf,[gName,'.png'],'png');
saveas(gcf,[gName,'.eps'],'psc2');
 
figure(9);  % fraction of covered, uncovered and overlap points regarding the size of the giant component 
 
plot(data(:,1),data(:,16),'LineStyle','-','LineWidth',2,'Marker','x');
hold all;
plot(data(:,1),data(:,17),'LineStyle','-','LineWidth',2,'Marker','o');
hold all;
plot(data(:,1),data(:,18),'LineStyle','-','LineWidth',2,'Marker','*');

title(strTitle,'Interpreter','latex','FontSize',14);

h = legend('f_{Spoints_c}','f_{Spoints_u}','Spoints_{overlap}','Interpreter','latex','Location','Best');

ylabel('Fraction f(S(g))/f(g)','FontSize',14);
xlabel(param.labelX,'FontSize',14);

annotation('textbox',[0.684928571428568 0.816666666666672 0.191071428571429 0.0690476190476191],...
    'String',strcat('Gains: \sigma=',sprintf('%d',param.gainConnectivityController),',\psi=',sprintf('%d',param.gainRobustnessControl),',\zeta=',sprintf('%d', param.gainCoverageController) ),...
    'FontSize',12,...
    'FontName','Cambria',...
    'FitBoxToText','off',...
    'LineStyle','none');


gName=[param.pathR,'graphic_Fraction_S_G_',param.fileName,'_',num2str(param.networkSize)];

saveas(gcf,[gName,'.fig'],'fig');
saveas(gcf,[gName,'.png'],'png');
saveas(gcf,[gName,'.eps'],'psc2');


close all;