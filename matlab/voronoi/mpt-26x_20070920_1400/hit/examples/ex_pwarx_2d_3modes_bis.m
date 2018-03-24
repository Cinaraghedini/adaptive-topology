%EX_PWARX_2D_3MODES_BIS Identification of a PWARX system with 3 modes from a raw dataset
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

% load the dataset (matrix Xid of regressors [y(k), u(k)]and vector yid of
% output samples y(k+1)

load data_PWARX_2d

% Setup fields of idpar. See hit_init for a description of the fields

% number of modes
idpar.s=3;

% size of Local Datasets
idpar.c=6;

% pattern_recognition algorithm
idpar.patt_rec_algo='mrlp';

% 
% Identification of the PWARX model
% 

% The structure idmodes describes the PWA map.
% Type 'help hit_regression' for a description of all the outputs of
% hit_regression

[idmodes,F,xi,LDs,inliers]=hit_regression(Xid,yid);

%
% Plot the identified PWARX model
% 

% plot the model in figure 2
hit_plot_idmodes(Xid, yid,inliers,idmodes,2);
hold on
title('Identified PWARX model');
ylabel('{u(k-1)}','FontSize',16)
xlabel('{y(k-1)}','FontSize',16)
zlabel('{y(k)}','FontSize',16)
set(gca,'FontSize',13)
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

% If the regressor set has dimension >1 and a pattern recognition algorithm
% different from MRLP has been used, test if the mode regions cover the
% regressor set or if some "holes" have been left.
% No hole has been left if H is an empty polytope. Otherwise it is a
% polytope array describing the holes.

% H=hit_holes(idmodes,idmodes.Regressor_set)

