function plotPropertiesDatabase(properties,param)

data = mean(properties);

datastd = var(properties);

figure1=figure;
 
axes1 = axes('Parent',figure1,'XTickLabel',{'k'}, 'XTick',[1], 'XMinorTick','on','TickDir','out');

box(axes1,'on');
hold(axes1,'all');

barwitherr(datastd(:,7), data(:,7),'BarWidth',1,'Parent',axes1);

title(['N - ' num2str(param.networkSize)  ', Area - ' num2str(param.area) , ' Range - ' num2str(param.range) ', Network - ' num2str(param.network)] );

gName=[param.path ,'graphicDegree'];

saveas(gcf,[gName,'.fig'],'fig');
saveas(gcf,[gName,'.jpg'],'jpg');
saveas(gcf,[gName,'.eps'],'psc2');

close all;

figure1=figure;
 
axes1 = axes('Parent',figure1,'XTickLabel',{'$\theta$','S'''}, 'XTick',[1 2], 'XMinorTick','on','TickDir','out');

box(axes1,'on');
hold(axes1,'all');

barwitherr(datastd(:,8:9), data(:,8:9),'BarWidth',1,'Parent',axes1);

title(['N - ' num2str(param.networkSize)  ', Area - ' num2str(param.area) , ' Range - ' num2str(param.range) ', Network - ' num2str(param.network)] );

gName=[param.path ,'graphicRobustness'];

saveas(gcf,[gName,'.fig'],'fig');
saveas(gcf,[gName,'.jpg'],'jpg');
saveas(gcf,[gName,'.eps'],'psc2');

close all;

figure1=figure;
 
axes1 = axes('Parent',figure1,'XTickLabel',{'eGlobal', 'eLocal','$\lambda$'}, 'XTick',[1 2 3], 'XMinorTick','on','TickDir','out');

box(axes1,'on');
hold(axes1,'all');

barwitherr([datastd(:,1) datastd(:,4) datastd(:,5)] , [data(:,1) data(:,4) data(:,5)],'BarWidth',1,'Parent',axes1);

title(['N - ' num2str(param.networkSize)  ', Area - ' num2str(param.area) , ' Range - ' num2str(param.range) ', Network - ' num2str(param.network)] );

gName=[param.path ,'graphicProperties'];

saveas(gcf,[gName,'.fig'],'fig');
saveas(gcf,[gName,'.jpg'],'jpg');
saveas(gcf,[gName,'.eps'],'psc');

