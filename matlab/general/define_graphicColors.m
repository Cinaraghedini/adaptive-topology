
function [dataPlot] = define_graphicColors(data)

pos=find(data(:,5)==1);

dataPlot=[];

for i=pos
  dataPlot=[dataPlot; data(i,:)];
end    
