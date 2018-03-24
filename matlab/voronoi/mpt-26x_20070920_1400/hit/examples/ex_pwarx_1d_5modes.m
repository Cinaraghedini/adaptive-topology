%EX_PWARX_1D_5MODES Identification of a PWARX system with 5 modes
%
% -------------------------------------------------------------------------
% DESCRIPTION
% -------------------------------------------------------------------------
%
% This is a batch file. 
% 

% initialize the HIT and MPT toolboxes
hit_init

% define idpar and plotpar as global to make them visible
% in the workspace
global idpar plotpar

% ----------- FIRST PART
% Generate input and output datapoints
% 

% N: number of generated datapoints
N=100;

% specify the mode PVs. Note that two of them have a continuous junction.
th_1 =[1 .2];
th_2 = [-1 2] ;
th_3 = [1 -1];
th_4 = [0 2];
th_5 = [2 -10];

Theta={th_1,th_2,th_3,th_4,th_5};

% Define the regressor set equal to the interval [-4,8]
Regressor_set=polytope([1 ;-1 ],[8;4]); 

% Regions: polyhedra array specifying the true mode regions.
% They are the intervals [-inf,-1], [-1,2], [2,4], [4,6], [6,inf]
Regions=[polytope([1],[-1]),polytope([1;-1],[2; 1]);polytope([-1;1],[-2;4]);polytope([-1;1],[-4;6]);polytope([-1],[-6])];

% standard deviation and variance of the noise corrupting
% the output measurements
sigma_sq=0.01;
sigma=sqrt(sigma_sq);

% create N samples of the scalar input 
% u(k) unifomly distributed in Regressor_set
V=extreme(Regressor_set);
u=hit_sample_intervals({{V,N}}, Regressor_set); 

% initialize the vector of output samples  with zeros.
% This means that, for computing y(2), the 'initial state' y(1)=0 is used
y=zeros(1,N+1);

% simulate the PWA system to obtain the output data
for k=1:N
    point=[u(k)];
    y(k+1)=hit_pwa(Theta,Regions,point)+sigma*randn(1); % add noise to the data
end

% Build the vectors of datapoints that will be used for identification
%
% Each mode of the PWARX system is assumed to be in the form
% y(k)=theta_i * [y(k-1) ... y(k-na) u'(k-1) ... u'(k-nb)]' + constant
%
% Xid(k,:) is the i-th regressor, i.e. 
% [y(k-1) ... y(k-na) u'(k-1) ... u'(k-nb)] for k>max([na,nb])
% yid(k) is the corresponding output measurements
na=0;
nb=1;
[Xid,yid]=hit_pwarx_format_data(u,y,na,nb);

% plot the true model and the data in figure 1
hit_pwa_plot1d(Theta,Regions,Regressor_set,1,'u(k-1)','y(k)');
hold on
plot(Xid(:,1),yid,'+','MarkerSize',8);
% save the axis
axis_saved=axis;
title('True model and datapoints')
hold off

% ----------- SECOND PART
% Setup all the parameters for hit_regression
% 

% Setup fields of idpar. See hit_init for a description of the fields

% number of modes
idpar.s=5;

% size of Local Datasets
idpar.c=8;

% use SVC as pattern_recognition algorithm
idpar.patt_rec_algo='svc';

% Define the regressor set
idpar.Regressor_set=Regressor_set;

% Define the regressor set *for simulation* equal to the interval [-10,10]
idpar.Regressor_set_sim=polytope([1 ;-1 ],[10;10]);


% ----------- THIRD PART
% Identification od the PWARX model
% 

% The structure idmodes describes the PWARX model.
% Type 'help hit_regression' for a description of all the outputs of
% hit_regression
[idmodes,F,xi,LDs,inliers]=hit_regression(Xid,yid);

% ----------- FOURTH PART
% Plot the identified PWARX model
% 

% plot the model in figure 2 on the regression set for simulation
flag_bigregionsim=1;
hit_plot_idmodes(Xid, yid,inliers,idmodes,2,flag_bigregionsim);
hold on
title('Identified PWARX model on the regression set for simulation');
xlabel('u(k-1)')
ylabel('y(k)')
% put the same axis as in the plot of the original model
% axis(axis_saved)
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
