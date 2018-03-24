% filename: coverageEvaluation.m
% Purpose:  computes the spots that are covered, uncovered and overlapped
% by the network nodes by generating N uniformly distributed points into a
% square area based on the node positions. The number os spots is given by param.nrPoints=50(the higher this value, the greater the accuracy).  
% Input: 
% - position - x and y coordinates for each network agent
% - param - parametrization struct
% - options - graph options (set the parametrization for graph properties
% computation
% Output: 
% coverageIdx array containing 
% - coverageIdx (1,1) fraction of covered spots 
% - coverageIdx (1,2) normalized fraction of covered spots
% - coverageIdx (1,3) fraction of uncovered spots 
% - coverageIdx (1,4) normalized fraction of uncovered spots
% - coverageIdx (1,5) overlap index
% - coverageIdx (1,6) fraction of covered spots regarding the giant component 
% - coverageIdx (1,7) normalized fraction of covered spots regarding the giant component
% - coverageIdx (1,8) fraction of uncovered spots regarding the giant component
% - coverageIdx (1,9) normalized fraction of uncovered spots regarding the giant component
% - coverageIdx (1,10) overlap index regarding the giant component 
% - coverageIdx (1,11) fraction of covered spots by the fraction of covered spots of the giant component  
% - coverageIdx (1,12) fraction of uncovered spots by the fraction of uncovered spots of the giant component
% - coverageIdx (1,13) fraction of overlapping spots by the fraction of overlapping spots of the giant component

function  [coverageIdx] =  coverageEvaluation(position,param,options)

XY = createPointsValidation(position,param); % generates spots into the environmet

A = computeDistancePoints(XY,position,param); % boolean matrix mapping nodes and spots with 1 when node covers that spot

[fractionInside, fractionOutside,fractionOverlap] = coverageEstimation(A); % computes the fraction of points XY covered, uncovered and overlapped by nodes

positionS = largestComponent(position,param,options); % find the largest connected component

XYS = createPointsValidation(positionS,param); % generates spots into the environmet

AS = computeDistancePoints(XYS,positionS,param); % boolean matrix mapping nodes and spots with 1 when node covers that spot
 
[fractionInsideS, fractionOutsideS, fractionOverlapS] = coverageEstimation(AS); % % computes the fraction of points XYS covered, uncovered and overlapped by nodes

normalizedFractionInside=fractionInsideS(1,1)/fractionInside(1,1);
normalizedFractionOutside=fractionOutsideS(1,1)/fractionOutside(1,1);
normalizedFractionOverlap=fractionOverlapS(1,1)/fractionOverlap(1,1);

coverageIdx=[fractionInside, fractionOutside,fractionOverlap,fractionInsideS, fractionOutsideS, fractionOverlapS,normalizedFractionInside,normalizedFractionOutside,normalizedFractionOverlap];