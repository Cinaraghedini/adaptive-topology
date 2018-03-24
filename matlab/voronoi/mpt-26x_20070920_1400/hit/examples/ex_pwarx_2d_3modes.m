%EX_PWARX_2D_3MODES Identification of a PWARX system with 3 modes
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
% Generate input and output datapoints
% 

% N: number of generated datapoints
N=60;

% specify the mode PVs. Note that two of them are equal
th_1 =[0.6 0.3 0];
th_2 = [-.6 -0.3 0];
th_3=[0.6 0.3 0];

% store the mode PVs in the cell array Theta
Theta={th_1,th_2,th_3};

% Define the regressor set equal to the square [-2,2]*[-2 2]
Regressor_set=polytope([1 0;-1 0;0 1; 0 -1],[2;2;2;2]); 

% Regions: polytope array specifying the true mode regions
Regions=[polytope([0.2 1],[-.8]),polytope([-0.2 1; -0.2 -1],[.8; .8]), polytope([0.2 -1],[-.8])];

% intersect Regions with the regressor set
for i=1:3
    Regions(i)=Regions(i)&Regressor_set;
end

% standard deviation and variance of the noise corrupting
% the output measurements
sigma_sq=0.01;
sigma=sqrt(sigma_sq);

% create N samples of the scalar input 
% u(k) unifomly distributed in [-2,2]
u=hit_sample_intervals({{[-2,2],N}}, polytope([1;-1],[2; 2])); 

% initialize the vector of output samples  
% using y(1)=0.5 as initial state

y=[0.5 zeros(1,N)]; 

% simulate the PWA system to obtain the output data
for k=1:N
    point=[y(k) u(k)]';
    [val,ind]=hit_pwa(Theta,Regions,point); 
    y(k+1)=val+sigma*randn(1); % add noise to the data
end

% Build regressors and outputs for the PWARX model of orders 1,1
% Xid(i,:) is the i-th regressors
% yid(i) is the i-th output sample
[Xid,yid]=hit_pwarx_format_data(u,y,1,1);

Nid=length(yid);

% x_grid, y_grid: grids for plotting the pwarx model
x_grid=-2:.1:2;
y_grid=-2:.1:2;

% plot the model and the data in figure 2
% plot the mode hyperplanes and the data in figure 2
minz=hit_pwa_plot2d(Theta,Regions,[],x_grid,y_grid,2,plotpar.color_surface);
hold on
for i=1:Nid
    hp(i) = plot3(Xid(i,1),Xid(i,2),yid(i),'ob');
    hold on
    set(hp(i), 'MarkerSize', 8,'LineWidth',1.5);
end
axis square
grid on
hold on
ylabel('{u(k-1)}','FontSize',16)
xlabel('{y(k-1)}','FontSize',16)
zlabel('{y(k)}','FontSize',16)
set(gca,'FontSize',13)
% plot the mode regions
hit_plot_regions3d(Regions,minz,plotpar.color_regions);
% save the axis
axis_saved=axis;
title('True PWA model');
hold off


% ----------- SECOND PART
% Set-up all the parameters for the identification
% algorithm.
% 

% Setup fields of idpar. See hit_init for a description of the fields

% number of modes
idpar.s=3;

% size of Local Datasets
idpar.c=6;

% regressor_set
idpar.Regressor_set=Regressor_set;

% pattern_recognition algorithm
idpar.patt_rec_algo='mrlp';

% ----------- THIRD PART
% The identification of the PWARX model
% 

% The structure idmodes describes the PWA map.
% Type 'help hit_regression' for a description of all the outputs of
% hit_regression

[idmodes,F,xi,LDs,inliers]=hit_regression(Xid,yid);

% ----------- FOURTH PART
% Plot the identified PWARX model
% 

% Setup the fields of plotpar
plotpar.x_grid=x_grid;
plotpar.x_grid=y_grid;

% plot the model in figure 3
hit_plot_idmodes(Xid, yid,inliers,idmodes,3);
hold on
% put the same axis as in the plot of the original model
axis(axis_saved)
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

% H=hit_holes(idmodes,idpar.Regressor_set)