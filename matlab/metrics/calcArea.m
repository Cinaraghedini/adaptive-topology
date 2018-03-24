%{
function calcArea: calculate coverage area of multirobots in a 2D space using geometric analysis 
    that builds inner and outer polygons from inner (blue) and outer (red)
    intersections of circunferences (range of each robot).
    The "outers" are external intersections and represented as red 'x' (when rendered by function plotScenario),
    each outer point closes a triangle with the centers of two envolved
    robots compounding a outer polygon.
    Inner polygons are composed by the bases of outer polygons in a region
    with inner intersections represented by blue '.'.
    Summing all outer polygons (red triangle), all inner (blue) polygons e
    all circular sector between each outer intersection leads to overall
    covered area.
    Run main.m file to see examples.

    Constraints: run 'init' function BEFORE to set path and related global variables.
    
input: robotsXYR: array having robots instances at lines, 
        x at column 1, 
        y at column 2 and
        circular range (radius) at column 3.

output: area with dimension of x-axis squared. 
    x-axis and y-axis must have same dimension.

@author Israel C. S. Rocha (israel.c.rocha@gmail.com)

%}

function [area, timeInSeconds] = calcArea(robotsXYR)
    
    global robots;
            
    robots = robotsXYR;
    tic
    evalIntersections();
    evalPolygons();
    area = calcPolygonsArea() + calcSectorsArea();
    timeInSeconds = toc;
    
end

    
    
