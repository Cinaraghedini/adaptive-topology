%EX_CAKE2 How to speed-up PWA regression with many datapoints.
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
% Generate datapoints
% 


% Nid: number of generated datapoints
Nid=1000;

% specify the mode PVs
th_1 =[4 2 3];
th_2 = [-6  6 -5];
th_3=[4 -2 -2];

% store the mode PVs in the cell array Theta
Theta={th_1,th_2,th_3};

% Define the regressor set equal to the square [-1,1]*[-1,1]
R=[1 0;-1 0;0 1; 0 -1];
r=[1;1;1;1];
Regressor_set=polytope(R,r); 

% Regions: polytope array specifying the true mode regions
Regions=[polytope([0 -1;-1 -1/sqrt(3);R],[0;0;r]),...
        polytope([1 1/sqrt(3);1 -1/sqrt(3);R],[0; 0;r]),...
        polytope([-1 1/sqrt(3);0 1;R],[0;0;r])];

% standard deviation and variance of the noise corrupting
% the output measurements
sig_sq=.1;
sigma=sqrt(sig_sq);

% Build the vectors of datapoints that will be used for identification
%
% Xid(i,:) is the i-th regressor
% yid(i) is the corresponding output measurements
const=0;
ampl=2;

Xid=ampl*(rand(Nid,2)-0.5)+const; 

yid=[];
for k=1:Nid
    point=[Xid(k,1) Xid(k,2)];
    [val,ind]=hit_pwa(Theta,Regions,point); % add noise to the data
    yid(k)=val+sigma*randn(1);
end
yid=yid(:);

% x_grid, y_grid: grids for plotting the pwa function
x_grid=-1:.05:1;
y_grid=-1:.05:1;

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
ylabel('{x_2}','FontSize',16)
xlabel('{x_1}','FontSize',16)
zlabel('{y}','FontSize',16)
set(gca,'FontSize',13)
% plot the mode regions
hit_plot_regions3d(Regions,minz,plotpar.color_regions);
% save the axis
axis_saved=axis;
title('True PWA model');
hold off

% ----------- SECOND PART
% Setup all the parameters for hit_regression
% 

% Setup fields of idpar. See hit_init for a description of the fields

% number of modes
idpar.s=3;

% size of Local Datasets
idpar.c=10;

% regressor_set
idpar.Regressor_set=Regressor_set;

% To speed up clustering with many datapoints, 
% lower the default number of Kmeans iterations (15) to 5
idpar.clustalgo.kmeans.repetitions=5; 

% To further speed up clustering you can also
% use scalar confidence measures by uncommenting the next 2 lines.
% However, clustering will be less precise.
%idpar.clustalgo.kmeans.init_centers='scalars'; 
%idpar.clustalgo.kmeans.centers='scalars';

% Do not plot classified xi-points during the running of hit_regression
plotpar.plot_class_xi_points_yn='N';

% To speed up the regression with many datapoints, 
% choose the fastest pattern_recognition algorithm
idpar.patt_rec_algo='psvc';

% ----------- THIRD PART
% PWA regression
% 

% The structure idmodes describes the PWA map.
% Type 'help hit_regression' for a description of all the outputs of
% hit_regression
[idmodes,F,xi,LDs,inliers]=hit_regression(Xid,yid);

% ----------- FOURTH PART
% If the regressor set is 1- or 2-dimensional
% the identified model can be plotted
% 

% plot the model in figure 2
hit_plot_idmodes(Xid, yid,inliers,idmodes,3);
hold on
% put the same axis as in the plot of the original model
axis(axis_saved)
title('Identified PWA model');
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
