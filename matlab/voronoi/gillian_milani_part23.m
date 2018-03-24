function [area thiessenWeigth thiessenAverage aritAverage]=gillian_milani_part23(x,y,bordery)

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
% or bordery={[0.01 8 10 10 7.01 2 0 0.01],[0 0 3 6 10 10 6 0]};


%Quick Initalization
%rainfall in [mm]

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
n=length(R);
[xInter yInter]=funny_intersections(bordery,s1,s2,m);

hold on
plot(xInter,yInter,'o','Color','red');
hold off

vect=linspace(0,1,n);
couleur=[1*vect' 0.5*vect'*1.6 1-1*vect'];

[R V]=intersection_allocation(V,R,xInter,yInter,s1,s2,bordery);

[R V]=bordery_vertices_allocation(V,R,bordery,s1,s2);

area1=area_V(V,R,couleur,bordery);

xdatacentre=s1{length(s1),:};
ydatacentre=s2{length(s2),:};
hold on
scatter(xdatacentre,ydatacentre,50,[0.5 0 0],'filled')
hold off


area=area1;
areaTot=sum(area1);
thiessenWeigth=area1/areaTot;
thiessenAverage=rainfall.*thiessenWeigth;
aritAverage=mean(rainfall);
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

function [output1 output2]=intersection_allocation(V,R,xInter,yInter,xdata,ydata,border)
xdatacentre=xdata{length(xdata),:};
ydatacentre=ydata{length(ydata),:};
n=length(xInter);
m=length(V);

for i=1:1:n %indice des intersections
    clear distance
    for j=1:1:m %indice des centres
        distance(j)=(xInter(i)-xdatacentre(j))^2+(yInter(i)-ydatacentre(j))^2;
    end
    distance=sort(distance);
    proche1(i)=distance(1);
    proche2(i)=distance(2);
    nvpoint(i,:)=[xInter(i) yInter(i)];
    indic(i)=m+i;
    V(indic,:)=nvpoint;
end

for j=1:1:m %indice des centres
   
    for i=1:1:n %indices des intersections
        IN = inpolygon(xInter(i),yInter(i),border{1},border{2});
        ON = inpolygon(xInter(i),yInter(i),border{1},border{2});
        if ON || ~IN
        distance=(xInter(i)-xdatacentre(j))^2+(yInter(i)-ydatacentre(j))^2;
        if sum(proche1(i)==distance)|| sum(proche2(i)==distance)
            R{j}=[R{j} indic(i)];
        else
        end
        end
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
v=length(xdata)-1;
lv=length(V);
m=length(xdatacentre);

for j=1:1:n-1 %indice des sommets de la bordure


    for k=1:1:m % indice de centre
        a={[]};
        for i=1:1:v %indice des lignes
            
            line1X=[X(j) xdatacentre(k)];
            line1Y=[Y(j) ydatacentre(k)];
            line2X=xdata{i};
            line2Y=ydata{i};
            a{i}=polyxpoly(line1X,line1Y,line2X,line2Y);

        end
        b=0;
        
        for i=1:1:v
            zt=a{i};
            if ~isempty(zt)
                b=1;
            end
        end
           if b==0
                numero(j)=k;
                hold on
                plot(X(j),Y(j),'X');
                hold off
            else
           end
    end
end

for j=1:1:n-1 %indice des sommets de la bordure
    
    indice=numero(j);
    w=lv+j;
    R{indice}=[R{indice} w];
    V(w,:)=[X(j) Y(j)];
end
end



function area1=area_V(V,R,couleur,border)
n=length(R);

for i=1:1:n
    test = sum((R{i}-1)~=0);
    
    if test > 2
        
            r=R{i};
            if r(1)==1
            r=r(2:length(r));
            else
                
            end
            xi=V(r,1);
            yi=V(r,2);
            
            IN = inpolygon(xi,yi,border{1},border{2});
            if min(IN)==1
            K = convhull(xi,yi);
            xi=xi(K);
            yi=yi(K);
            area1(i)=polyarea(xi,yi);
            hold on
            fill(xi,yi,couleur(i,:));
            hold off
            else
                xi1=1./IN.*xi(:);
                xi2=xi;
                yi2=yi;
                xi=[];
                yi=[];
                for tv=1:1:length(xi2)
                   if xi1(tv)~=inf
                     
                       xi=[xi xi2(tv)];
                       yi=[yi yi2(tv)];
                   end
                end
            K = convhull(xi,yi);
            xi=xi(K);
            yi=yi(K);
            area1(i)=polyarea(xi,yi);
            hold on
            fill(xi,yi,couleur(i,:));
            hold off
            end
    else
    end
end
end 