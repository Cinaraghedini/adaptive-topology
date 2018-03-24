function [area thiessenWeigth thiessenAverage]=gillian_milani(x,y,bordery,rainfall)

% GILLIAN_MILANI: draw the thiessen (or voronoi) polygons inside a border
%
%Outputs :[area thiessenWeigth thiessenAverage]
%'area' gives back a vector containing the area of each polygons ordrered
%as x
%'thiessenWeigth gives back a vector with the thiessen Weigths
%thiessenAverage gives back a vector with the rainfall estimated on each
%region, Enjoy!
%
%the x and y inputs are the corrdinates of points
%the 3 argument, the border, must be given in the form :
%border={[0 8 10 10 10 7.01 2 0 0],[0 0 3 6 10 10 10 6 0]} for example for
%a octogone.
%
%ATTENTION: please if you have a border vertex exactly on a
%line of polygon going to infinite : change the border vertex with a small
%factor (so the vertex don't stand on the line) Maybe this will be better
%on 2nd version ;-)
%
%ATTENTION: all center points must be given inside the border
%
%ATTENTION: the last point of the border must be same as the first
%
%ATTENTION: this function give FALSE result when some geometric
%complexities implies occur about the vertices of the border (see
%bordery_vertices_allocation nested function). Howeverm this can be checked
%on the graph.
%
%example of data :
% x = [2 8 5 5 6 9]';
% y = [5 4 5 9 1 9]';
% rainfall = [10 7 5 12 2 8];
% bordery={[0 0 10 10 0],[0 10 10 0 0]};
% 
% clear all
% clc
% close all

%creation of Thiessen polygons
length(x);
voronoi(x, y);
dt = DelaunayTri([x y]);
[V,R] = voronoiDiagram(dt);

h = findobj(gca,'Type','line');
set(h,'Color','k')

%plotting of the boundary
xAxMin=min(bordery{1})-2;
xAxMax=max(bordery{1})+2;
yAxMin=min(bordery{2})-2;
yAxMax=max(bordery{2})+2;

hold on
plot(bordery{1},bordery{2});
hold off
axis([xAxMin xAxMax yAxMin yAxMax]);

%Master JOB!!! lokking for the intersection...
%between lines and boundary


s1=get(h,'XData');
s2=get(h,'YData');

m=length(s1)-1;
[xInter yInter]=funny_intersections(bordery,s1,s2,m);

hold on
plot(xInter,yInter,'o','Color','red');
hold off

n=length(R);
vect=linspace(0,1,n);
couleur=[1*vect' 0.5*vect'*1.6 1-1*vect'];

[R V]=intersection_allocation(V,R,xInter,yInter);
[R V]=bordery_vertices_allocation(V,R,bordery,s1,s2);

area1=area_V(V,R,couleur);

% newX=s1(1:m)
% newY=s2(1:m);

% set(findobj(gca,'Type','line'),'XData',newX)
% set(findobj(gca,'Type','line'),'YData',newY)


% h = findobj(gca,'Type','line');
% s1=get(h,'XData');
% s2=get(h,'YData');
%
% result=[s1 s2];
% area=area_polygon(V, R, axisVertices);
area=area1;
areaTot=sum(area1);
thiessenWeigth=area1/areaTot;
thiessenAverage=rainfall.*thiessenWeigth;
end

function [newX newY]=funny_intersections(bord,xdata,ydata,m)
j=1;
for i=1:1:m
    x=xdata{i};
    y=ydata{i};
    [xNew yNew]=polyxpoly(x, y, bord{1}, bord{2});
    if isempty(xNew)||isempty(yNew)
    else
    newX(j)=xNew;
    newY(j)=yNew;
    j=j+1;
    end
end
end

function [output1 output2]=intersection_allocation(V,R,xInter,yInter)

n=length(xInter);
m=length(V);
for i=1:1:n
    clear distance
    for j=1:1:m
        
        distance(j)=(xInter(i)-V(j,1))^2+(yInter(i)-V(j,2))^2;
        if distance(j)==min(distance);
            proche(i)=j;
            proche2(i,:)=[xInter(i) yInter(i)];
            proche3(i)=m+i;
            V(proche3,:)=proche2;
        end
    end
end


for i=1:1:m %indice des polygones
    y=R{i};
    if y(1)==1
    for j=1:1:n %indices des intersections
        if sum(proche(j)==y)~=0
            R{i}=[R{i} proche3(j)];
        else
        end
    end
      
    else
    end
end
output1=R;
output2=V;
end


function [R V] = bordery_vertices_allocation(V,R,bordery,xdata,ydata)
%but de la fonction : inclure les sommet de la bordure dans les polygones
%respectif
X=bordery{1};
Y=bordery{2};
xdatacentre=xdata{length(xdata),:};
ydatacentre=ydata{length(ydata),:};
n=length(X);
v=length(V);
m=length(xdatacentre);
for j=1:1:n-1 %indice des sommets de la bordure
    %cette boucle est sensé assigné le centre de polygone le plus proche
    %a chaque sommet de la bordure
    z=0;
    for i=1:1:m %indice de centre
        %cette boucle est sensé calculé le centre le plus proche au sommet X(j)
        %Y(j)
        z(i) = (xdatacentre(i)-X(j))^2+(ydatacentre(i)-Y(j))^2;
        hold on
        plot(X(j),Y(j),'X')
        hold off
        if z(i) == min(z);
           numero(j)=i;
        else
        end
    end
end
for j=1:1:n-1 %indice des sommets de la bordure
    indice=numero(j);
    w=v+j;
    R{indice}=[R{indice} w];
    V(w,:)=[X(j) Y(j)];
end

end



function area1=area_V(V,R,couleur)
n=length(R);

for i=1:1:n
    test = sum((R{i}-1)~=0);
    
    if test > 2
        y=R{i};
        
        if y(1)==1
            r=R{i};
            r=r(2:length(r));
            xi=V(r,1);
            yi=V(r,2);
            K = convhull(xi,yi);
            xi=xi(K);
            yi=yi(K);
            area1(i)=polyarea(xi,yi);
            hold on
            fill(xi,yi,couleur(i,:));
            hold off
        else
            xi=V(R{i},1);
            yi=V(R{i},2);
            area1(i)=polyarea(xi,yi);
            hold on
            fill(xi,yi,couleur(i,:));
            hold off
        end
        
    else
    end
end
end

