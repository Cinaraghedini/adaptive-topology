% filename: compute_a.m
% Purpose: compute a coefficient for compousing a quadratic equation as defined in the referentece
% Input:
% - v1 - array with x and y coordinates of vertice v1 of a line
% - v2 - array with x and y coordinates of vertice v2 of a line
% Output:
% a - coefficient a value
% - Reference :Erik Tjong Kim Sang. Voronoi diagrams without bounding boxes. In ISPRS Annals of the Photogrammetry, 
% Remote Sensing and Spatial Information Sciences, II-2/W2, 2015.

function [a] = compute_a(v1,v2)
if v1(1,1)==v2(1,1)
    a=0;
else
    a=(v2(1,2)-v1(1,2))/(v2(1,1)-v1(1,1));
end