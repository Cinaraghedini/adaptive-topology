% filename: compute_b.m
% Purpose: compute b coefficient for compousing a quadratic equation as defined in the referentece
% Input:
% - v1 - array with x and y coordinates of vertice v1 of a line
% - v2 - array with x and y coordinates of vertice v2 of a line
% Output:
% b - coefficient b value
% - Reference :Erik Tjong Kim Sang. Voronoi diagrams without bounding boxes. In ISPRS Annals of the Photogrammetry, 
% Remote Sensing and Spatial Information Sciences, II-2/W2, 2015.

function [b] = compute_b(v1,v2)

y1=v1(1,2);
y2=v2(1,2);

x1=v1(1,1);
x2=v2(1,1);

if x1==x2
    b=x1;
else
    b=(y1*x2-y2*x1)/(x2-x1);
end

