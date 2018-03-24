%EX_APPROX1D Approximate a sin function with a continuous PWA function
%
% -------------------------------------------------------------------------
% DESCRIPTION
% -------------------------------------------------------------------------
%
% This is a batch file.
% 

% Copyright is with the following author:
%
% (C) 2005 Giancarlo Ferrari Trecate,
%         giancarlo.ferrari@unipv.it
% -------------------------------------------------------------------------


% initialize the HIT and MPT toolboxes
hit_init

% define idpar and plotpar as global to make them visible
% in the workspace
global idpar plotpar

% ----------- FIRST PART
% Sample the function

% N: number of generated datapoints
N=60 ;

% create a grid in the function domain
Xid=[-pi:2*pi/N:pi];

% function samples
yid=sin(Xid);

% Regressor_set: interval [-pi,pi]
idpar.Regressor_set=polytope([1 ;-1 ],[pi;pi]); 

% plot the data in figure 1

figure(1);clf;
hold on
plot(Xid,yid,'+','MarkerSize',8);
hold on
baseline=axis;
title('Function samples')
hold off

% ----------- SECOND PART
% Setup all the parameters for hit_regression
% 

% Setup fields of idpar. See hit_init for a description of the fields

% Number of modes of the PWA approximator
idpar.s=5;
% Size of Local Datasets
idpar.c=8;

% pattern_recognition algorithm
idpar.patt_rec_algo='svc';

% we want a continuous PWA approximation
idpar.continuity='c';

% ----------- THIRD PART
% PWA regression
% 

% The structure idmodes describes the PWA map.
% Type 'help hit_regression' for a description of all the outputs of
% hit_regression

% put Xid in column format
Xid=Xid(:);
[idmodes,F,xi,LDs,inliers]=hit_regression(Xid,yid);

% ----------- FOURTH PART
% If the regressor set is 1- or 2-dimensional
% the PWA model can be plotted
% 

% avoid plotting datapoints
plotpar.datapoints_yn='N';
hit_plot_idmodes(Xid, yid,inliers,idmodes,2);

% plot the original function on the same figure
figure(2);
hold on
plot(Xid,yid,'r');
title('True function and PWA approximation');
hold off

% ----------- LAST PART
% validate the model
% 

% Test if all modes have a "similar" mse.
% Modes with "high" mse are likely to be badly reconstructed. 

% [mse,mse_mode]=hit_mse(Xid,yid,idmodes);

% Test if the pattern recognition algorithm has misclassified few
% datapoints.
% Type 'help hit_regression' for the meaning of 
% idmodes.pattern_rec_valid.correctness
% A correctness of 100 means that no point has been misclassified.

% idmodes.pattern_rec_valid.correctness

% If you used Kmeans for clustering, check if Kmeans converged to the
% minimal cost in many different runs. Otherwise, either clusters are not
% well-separated in the xi-space (in this case increase idpar.c) or Kmeans
% started from a bad random initialization in all runs (in this case
% increase idpar.clustalgo.kmeans.repetitions).
% Look into hit_init.m for the meaning of 
% idpar.clustalgo.kmeans.repetitions

% plot(idmodes.clust_valid.costs,'*')