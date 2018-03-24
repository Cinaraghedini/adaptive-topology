%EX_ESTIMATE_C_S Cross-validation for estimating the mode number and the LD size
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
% Generate regression and validation datapoints
% 

% N: number of datapoints
N=150 ;
% Nv: number of validation datapoints
Nv=N/2 ;

% specify the mode PVs: 
th_1 =[1 .2];
th_2 = [-1 2] ;
th_3 = [1 .2];
th_4 = [3 -10];
th_5 = [-2 15];

% store the mode PVs in the cell array Theta
Theta={th_1,th_2,th_3,th_4,th_5};

% Regions: polyhedra array specifying the true mode regions.
% They are the intervals [-inf,-1], [-1,2], [2,4], [4,6], [6,inf]
Regions=[polytope([1],[-1]),polytope([1;-1],[2; 1]);polytope([-1;1],[-2;4]);polytope([-1;1],[-4;6]);polytope([-1],[-6])];

% Define the regressor set equal to the interval [-4,8]
idpar.Regressor_set=polytope([1 ;-1 ],[8;4]);

% standard deviation and variance of the noise corrupting
% the output measurements
sigma_sq=0.1;
sigma=sqrt(sigma_sq);


% create input samples
% u(k) unifomly distributed in Regressor_set
V=extreme(idpar.Regressor_set);
u=hit_sample_intervals({{V,N}}, idpar.Regressor_set);

% create validation input
uv=hit_sample_intervals({{V,Nv}}, idpar.Regressor_set);

% create the output samples
y=zeros(1,N);
yv=zeros(1,Nv);

% sample the PWA map to obtain output samples
for k=1:N
    point=[u(k)];
    y(k)=hit_pwa(Theta,Regions,point)+sigma*randn(1);
end

% sample the PWA map to obtain output validation samples
for k=1:Nv
    point=[uv(k)];
    yv1(k)=hit_pwa(Theta,Regions,point)+sigma*randn(1);
end

% Identification regressors and outputs 
Xid=u(:);
yid=y(:);

% Validation regressors
Xv=uv(:);
yv=yv1(:);

% plot the true model and the data in figure 1
hit_pwa_plot1d(Theta,Regions,idpar.Regressor_set,1,'x(k)','y(k)');
hold on
plot(Xid(:,1),yid,'+','MarkerSize',8);
hold on
baseline=axis;
title('True model and datapoints')
hold off

% plot the true model and the validation data in figure 2
hit_pwa_plot1d(Theta,Regions,idpar.Regressor_set,2,'x(k)','y(k)');
hold on
plot(Xv(:,1),yv,'+','MarkerSize',8);
hold on
baseline=axis;
title('True model and validation datapoints')
hold off

% ----------- THIRD PART
% Estimate the best Local Dataset size and the best number of modes  
c_test=[6 8];
s_test=[4 5 6];

% To speed up the regressions made in hit_estimate_cs
% choose the fastest pattern_recognition algorithm
idpar.patt_rec_algo='psvc';
% Do not plot classified xi-points during the running of hit_estimate_cs
plotpar.plot_class_xi_points_yn='N';
[best_s,best_c,mse_m]=hit_estimate_cs(Xid,yid,Xv,yv,c_test,s_test);

% ----------- FOURTH PART
% PWA regression with the best Local Dataset size and the best number of
% modes 
idpar.c=best_c;
idpar.s=best_s;
% Plot classified xi-points during the running of hit_regression
plotpar.plot_class_xi_points_yn='y';
% Choose a more precise pattern_recognition algorithm for getting the final
% model
idpar.patt_rec_algo='svc';
[idmodes,F,xi,LDs,inliers]=hit_regression(Xid,yid);

% ----------- FIFTH PART
% Plot the identified model in figure 3
hit_plot_idmodes(Xid, yid,inliers,idmodes,3);
xlabel('x(k)')
ylabel('y(k)')

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
