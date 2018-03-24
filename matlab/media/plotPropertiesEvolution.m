% filename: plotPropertiesEvolution.m
% Purpose:  plot and persist network properties evolution for a specific
% gain configuration
% Input: 
% - data - matrix containing the coverage data evoluation
% - state - array with the average disconnection time
% - param - parametrization struct


function plotPropertiesEvolution(data,state,param)

close all;

strTitle=sprintf('N=%d, Area=$%d^2$, Iterations=%d', param.networkSize,param.area,param.numberNetworks);

figure(1);  % global efficiency
 
plot(data(:,1),data(:,2),'LineStyle','-','LineWidth',2,'Marker','x');

title(strTitle,'Interpreter','latex','FontSize',14);
ylabel('Global efficiency','FontSize',14);
xlabel('t','FontSize',14);
ylim([0 1]); 

annotation('textbox',[0.684928571428568 0.816666666666672 0.191071428571429 0.0690476190476191],...
    'String',strcat('Gains: \sigma=',sprintf('%d',param.gainConnectivityController),',\psi=',sprintf('%d',param.gainRobustnessControl),',\zeta=',sprintf('%d', param.gainCoverageController) ),...
    'FontSize',12,...
    'FontName','Cambria',...
    'FitBoxToText','off',...
    'LineStyle','none');


gName=[param.pathR,'graphic_eG_',param.fileName,'_',num2str(param.networkSize)];

saveas(gcf,[gName,'.fig'],'fig');
saveas(gcf,[gName,'.png'],'png');
saveas(gcf,[gName,'.eps'],'psc2');

figure(2);  % local efficiency
 
plot(data(:,1),data(:,3:5),'LineStyle','-','LineWidth',2,'Marker','x');

title(strTitle,'Interpreter','latex','FontSize',14);

ylabel('Local efficiency','FontSize',14);
xlabel('t','FontSize',14);
ylim([0 1]); 

h = legend('Local Efficiency','Clustering Coefficient','Weighted Cluter Coefficient','Location','Best');

annotation('textbox',[0.684928571428568 0.816666666666672 0.191071428571429 0.0690476190476191],...
    'String',strcat('Gains: \sigma=',sprintf('%d',param.gainConnectivityController),',\psi=',sprintf('%d',param.gainRobustnessControl),',\zeta=',sprintf('%d', param.gainCoverageController) ),...
    'FontSize',12,...
    'FontName','Cambria',...
    'FitBoxToText','off',...
    'LineStyle','none');


gName=[param.pathR,'graphic_eL_', param.fileName,'_', num2str(param.networkSize)];

saveas(gcf,[gName,'.fig'],'fig');
saveas(gcf,[gName,'.png'],'png');
saveas(gcf,[gName,'.eps'],'psc2');


figure(3);  % giant component
 
plot(data(:,1),data(:,7),'LineStyle','-','LineWidth',2,'Marker','x');

title(strTitle,'Interpreter','latex','FontSize',14);

ylabel('S','FontSize',14);
xlabel('t','FontSize',14);
ylim([0 1]); 

annotation('textbox',[0.684928571428568 0.816666666666672 0.191071428571429 0.0690476190476191],...
    'String',strcat('Gains: \sigma=',sprintf('%d',param.gainConnectivityController),',\psi=',sprintf('%d',param.gainRobustnessControl),',\zeta=',sprintf('%d', param.gainCoverageController) ),...
    'FontSize',12,...
    'FontName','Cambria',...
    'FitBoxToText','off',...
    'LineStyle','none');

gName=[param.pathR,'graphic_S_',param.fileName ,'_',num2str(param.networkSize)];

saveas(gcf,[gName,'.fig'],'fig');
saveas(gcf,[gName,'.png'],'png');
saveas(gcf,[gName,'.eps'],'psc2');
    

figure(4);  % lambda - algebric connectivity
 
plot(data(:,1),data(:,6),'LineStyle','-','LineWidth',2,'Marker','x');

title(['N - ' num2str(param.networkSize)  ', Area - ' num2str(param.area) , ', Iterations - ' num2str(param.numberNetworks)] ,'FontSize',14);

ylabel('\lambda','FontSize',14);
xlabel('t','FontSize',14);

annotation('textbox',[0.684928571428568 0.816666666666672 0.191071428571429 0.0690476190476191],...
    'String',strcat('Gains: \sigma=',sprintf('%d',param.gainConnectivityController),',\psi=',sprintf('%d',param.gainRobustnessControl),',\zeta=',sprintf('%d', param.gainCoverageController) ),...
    'FontSize',12,...
    'FontName','Cambria',...
    'FitBoxToText','off',...
    'LineStyle','none');

gName=[param.pathR,'graphic_Ag_',param.fileName,'_', num2str(param.networkSize)];

saveas(gcf,[gName,'.fig'],'fig');
saveas(gcf,[gName,'.png'],'png');
saveas(gcf,[gName,'.eps'],'psc2');
    
figure(5);  % Average Degree - k 
 
plot(data(:,1),data(:,8),'LineStyle','-','LineWidth',2,'Marker','x');

title(strTitle,'Interpreter','latex','FontSize',14);

ylabel('<k>','FontSize',14);
xlabel('t','FontSize',14);
annotation('textbox',[0.684928571428568 0.816666666666672 0.191071428571429 0.0690476190476191],...
    'String',strcat('Gains: \sigma=',sprintf('%d',param.gainConnectivityController),',\psi=',sprintf('%d',param.gainRobustnessControl),',\zeta=',sprintf('%d', param.gainCoverageController) ),...
    'FontSize',12,...
    'FontName','Cambria',...
    'FitBoxToText','off',...
    'LineStyle','none');


gName=[param.pathR,'graphic_K_',param.fileName,'_',num2str(param.networkSize)];

saveas(gcf,[gName,'.fig'],'fig');
saveas(gcf,[gName,'.png'],'png');
saveas(gcf,[gName,'.eps'],'psc2');


figure(6);  % Robustness level
 
plot(data(:,1),data(:,9),'LineStyle','-','LineWidth',2,'Marker','x');

title(strTitle,'Interpreter','latex','FontSize',14);

ylabel('\Theta','FontSize',14);
xlabel('t','FontSize',14);

ylim([0 1]); 

annotation('textbox',[0.684928571428568 0.816666666666672 0.191071428571429 0.0690476190476191],...
    'String',strcat('Gains: \sigma=',sprintf('%d',param.gainConnectivityController),',\psi=',sprintf('%d',param.gainRobustnessControl),',\zeta=',sprintf('%d', param.gainCoverageController) ),...
    'FontSize',12,...
    'FontName','Cambria',...
    'FitBoxToText','off',...
    'LineStyle','none');

gName=[param.pathR,'graphic_R_',param.fileName,'_',num2str(param.networkSize)];

saveas(gcf,[gName,'.fig'],'fig');
saveas(gcf,[gName,'.png'],'png');
saveas(gcf,[gName,'.eps'],'psc2');


figure(7);  % giant component generated by network diconnection during the robustness evaluation 
 
plot(data(:,1),data(:,10),'LineStyle','-','LineWidth',2,'Marker','x');

title(strTitle,'Interpreter','latex','FontSize',14);

ylabel('S - \Theta','FontSize',14);
xlabel('t','FontSize',14);
ylim([0 1]); 

annotation('textbox',[0.684928571428568 0.816666666666672 0.191071428571429 0.0690476190476191],...
    'String',strcat('Gains: \sigma=',sprintf('%d',param.gainConnectivityController),',\psi=',sprintf('%d',param.gainRobustnessControl),',\zeta=',sprintf('%d', param.gainCoverageController) ),...
    'FontSize',12,...
    'FontName','Cambria',...
    'FitBoxToText','off',...
    'LineStyle','none');

gName=[param.pathR,'graphic_SR_',param.fileName,'_',num2str(param.networkSize)];

saveas(gcf,[gName,'.fig'],'fig');
saveas(gcf,[gName,'.png'],'png');
saveas(gcf,[gName,'.eps'],'psc2');

if ~ isempty(state) % average disconnection time
    figure(8);
    
    bar(state(1,:),'DisplayName','dds(:,1)');
    
    title(['N - ' num2str(param.networkSize)  ', Area - ' num2str(param.area) , ', Iterations - ' num2str(param.numberNetworks)] ,'FontSize',14);
    
    xlabel('Iteration time','FontSize',14);
    ylim([0 param.tf]);
    
    gName=[param.pathR,'graphic_ItT_',param.fileName,'_',num2str(param.networkSize)];
    
    saveas(gcf,[gName,'.fig'],'fig');
    saveas(gcf,[gName,'.jpg'],'jpg');
    saveas(gcf,[gName,'.eps'],'psc2');
end
close all;