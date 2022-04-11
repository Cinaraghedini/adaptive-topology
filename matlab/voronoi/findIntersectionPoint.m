% filename: findIntersectionPoint.m
% Purpose:  finding the intersection point given a vertice and a set of
% candidate points
% Input: 
% - seedPoints : matrix of size N containing the candidates points for being
% the intersection point
% - posPoint: the position of the target point 
% Output:
% - xy  - x and y coordinates of the intersection point 
% - Reference :Erik Tjong Kim Sang. Voronoi diagrams without bounding boxes. In ISPRS Annals of the Photogrammetry, 
% Remote Sensing and Spatial Information Sciences, II-2/W2, 2015.

function [xy] = findIntersectionPoint(seedPoints,posPoint)

for i=1:size(seedPoints,1)

    x=[seedPoints(i,1) posPoint(1,1)];
    y=[seedPoints(i,2) posPoint(1,2)];
     
    % finds the coefficients of a polynomial P(x) of degree 1 that fits the data 1 best
    fx(i,:)=polyfit(x,y,1); 
    
    % compute the distance between the point and the possible points
%    try
      d(i,:)=pdist(real([seedPoints(i,:); posPoint]),'euclidean');
  %      d(i,:)=pdist([seedPoints(i,:); posPoint],'euclidean');
    
       % disp(seedPoints(i,:))
       % disp(posPoint)
 %   catch
 %       disp('a')
 %    
 %       disp(real([seedPoints(i,:); posPoint]))
 %       disp([seedPoints(i,:); posPoint])
        
    
        % disp(seedPoints(i,:))
       % disp(posPoint)
  %  end
end

% the coefficient are equal for the roots, defines the closest one as
% intersection point, otherwise the higher coeficient
if abs(fx(1,1)-fx(2,1)) < 10^-3
    pos=find(d(:,1)==min(d(:,1)));
else
    pos=find(fx(:,1)==max(fx(:,1)));
end
xy=seedPoints(pos,:);
