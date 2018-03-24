function hit_init
%HIT_INIT Initializes the HIT and MPT toolboxes
%
% -------------------------------------------------------------------------
% DESCRIPTION
% -------------------------------------------------------------------------
% 
% Many routines of the HIT and MPT toolboxs can be called with 
% user-specified values of different parameters. 
% Default values of the HIT toolbox are stored in the variables idpar and
% plotpar, 
% Default values of the MPT toolbox are stored in the variable mptOptions
% All of them are are stored in MATLAB's workspace as a global variables. 
% 
% Please have a look at the code of hit_init.m to see what are the default 
% values and the meaning of the fields
%
% -------------------------------------------------------------------------
% INPUT
% -------------------------------------------------------------------------
% none
%
% -------------------------------------------------------------------------
% OUTPUT                                                                                                    
% -------------------------------------------------------------------------
% global variables idpar, plotpar and mptOptions
%

% Copyright is with the following author:
%
% (C) 2005 Giancarlo Ferrari Trecate,
%         giancarlo.ferrari@unipv.it
% -------------------------------------------------------------------------
% Legal note:
%     This program is free software; you can redistribute it and/or
%     modify it under the terms of the GNU General Public
%     License as published by the Free Software Foundation; either
%     version 2.1 of the License, or (at your option) any later version.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%     General Public License for more details.
%
%     You should have received a copy of the GNU General Public
%     License along with this library; if not, write to the
%     Free Software Foundation, Inc.,
%     59 Temple Place, Suite 330,
%     Boston, MA  02111-1307  USA
%
%
% -------------------------------------------------------------------------

% initialize the MPT toolbox
global mptOptions idpar plotpar;

if ~isstruct(mptOptions)
    mpt_init;
end


% Set the toolbox version

hit_ver='1.00';

% NEXT, THE STRUCTURE IDPAR IS DEFINED.
% CHANGE ONLY THE VALUES BUT *NOT* THE NAME OF THE FIELDS.

% Number of submodels. The empty matrix means that a value *must* be
% supplied by the user, if a supervised clustering algorithm 
% (like Kmeans) will be used

idpar.s=[];

% Size of Local Datasets (LDs). The empty matrix means that a value *must*
% be supplied by the user

idpar.c=[];

% idpar.Regressor_set: polytope specifying the regressor set 
%   (the big set that contains the regions). If it is empty and it is not
%   supplied by the user before calling hit_regression, it is set
%   automatically as the smallest hypercube containing all regressors.

idpar.Regressor_set=[]; 

% idpar.Regressor_set_sim: polytope  specifying the regressor set 
%   used for simulating the model. It is optional.
%   In some experiments one may want to identify a model on a "small"
%   Regressor_set and simulate the model over a bigger Regressor_set_sim.
%   If it is empty and it is not
%   supplied by the user before calling hit_regression, it is understood
%   that it will be not used.

idpar.Regressor_set_sim=[];

% idpar.Weight_primal: vector with scalar weights for datapoints. If
% nonepty (or nonexistent) data in LDs are weighted before computing
% Local Parameter Vectors (LPVs). This is an experimental feature.
% Default: []

idpar.Weight_primal=[];

% idpar.patt_rec_algo: 'MRLP'; 'SVC'; 'PSVC' 
% Specify the pattern recognition algorithm for finding the regions
% 'MRLP' = Multicategory Robust Linear Programming (LP-based. It gives the
% best results but it can be very slow)
% 'SVC' = Support Vector Classification (QP-based, pairwise separation)
% 'PSVC' = Proximal Support Vector Classification (Pairwise
% separation. It is very fast but it is also the less accurate algorithm.)
%
% Please note that if the regressor set has more than one dimension, SVC
% and PSVC may leave 'holes' in the regressor set, while MRLP partitions
% the whole regressor set.
% 
% Default: 'PSVC'

idpar.patt_rec_algo='PSVC';

% Choosing the solver for Linear Programming (LP) problems
% idpar.LPsolver, and idpar.LPsolver_cnstr_reduction can take one of the
% following values
%                0: uses E04MBF.M
%                1: uses linprog.m
%                2: uses CPLEX
%                3: uses CDD Criss-Cross
%                4: uses GLPK
%                5: uses CDD Dual Simplex
%                6: uses SeDuMi
%                7: uses QSopt
%
% LP problems are solved through the function mpt_solvelp of the MPT
% toolbox. See the manual of the MPT toolbox for various available solvers
%
% idpar.LPsolver specifies the solver to be used for linear programming in
% MRLP (default: mptOptions.lpsolver)
% idpar.LPsolver_cnstr_reduction specifies the solver to be used for
% removing redundant constraints in polytopes 
%
% Default:mptOptions.lpsolver

idpar.LPsolver=mptOptions.lpsolver; 
idpar.LPsolver_cnstr_reduction=mptOptions.lpsolver; 

% Choosing the solver for Quadratic Programming (QP) problems
% idpar.QPsolver specifies the solver to be used for finding models with
% the continuity constraints or for using SVC.
% It can take one of the following values
%                0:  uses E04NAF.M
%                1:  uses QUADPROG.M
%                2:  uses CPLEX
% QP problems are solved through the function mpt_solveqp of the MPT
% toolbox. See the manual of the MPT toolbox for various available solvers
%
% Default: mptOptions.qpsolver

idpar.QPsolver=mptOptions.qpsolver; 

% idpar.continuity: 'C', 'D'
% specify if the identified model must be continuous ('C')
% or discontinuous ('D')
%
% Default: 'D'

idpar.continuity='D'; 

% idpar.what_to_clust: 'LPVS' | 'FVS'
% Specify if Local Parameter Vectors (LPVs) or Feature Vectors (FVs) must
% be clustered.
%
% Default:FVs

idpar.what_to_clust='FVs';

% idpar.clustalgo.name: 'KMEANS' | 'SL'
% Specify the clustering algorithm. 
% 'KMEANS' is a modified KMEANS
% 'SL' is the single-linkage procedure discussed in
%       G. Ferrari-Trecate and M. Muselli. Single-linkage clustering 
%       for optimal classification in piecewise affine regression. 
%       IFAC Conference on the Analysis and Design of Hybrid Systems (ADHS 03), 2003.
%       It estimates automatically the number of modes

idpar.clustalgo.name='KMEANS'; 

% idpar.patt_rec.remove_duplicates: 'Y', 'N'
% Specify if LPVs/FVs duplicated (i.e. generated from
% identical LDs) must be removed before clustering
% 

idpar.clustalgo.remove_duplicates='N'; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETERS FOR Kmeans
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% n. of times Kmeans is repeated with different
% initializations.
idpar.clustalgo.kmeans.repetitions=15; 

% idpar.clustalgo.kmeans.init_centers: 'COVARIANCES'; 'SCALARS'
% Specify how to weight the xi-points for initializing 
% the centers of the clusters.
%'SCALARS': Use the peaks of Gaussian associated to each xi-point as
% weight.
%'COVARIANCES': Use the covariances associated to each xi-point as weight.
%
% Default: COVARIANCES

idpar.clustalgo.kmeans.init_centers='COVARIANCES'; 

% idpar.clustalgo.centers:  the same as idpar.clustalgo.init_centers
% but here one chooses the weights used for updating the 
% centers in Kmeans
%
% Default: COVARIANCES
idpar.clustalgo.kmeans.centers='COVARIANCES'; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END PARAMETERS FOR Kmeans
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETERS FOR Single-Linkage (SL)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% idpar.clustalgo.sl.guess_min_dist: minimal distance between clusters
% Used in SL clustering

idpar.clustalgo.sl.guess_min_dist=[];

% idpar.clustalgo.sl.plot_steps: 'Y','N'. If 'Y', each step of SL is plot.
% in the figure number idpar.clustalgo.sl.plot_steps.fig

idpar.clustalgo.sl.plot_steps='N';
idpar.clustalgo.sl.plot_fig=12;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END PARAMETERS FOR Single-Linkage (SL)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% idpar.mix_detect:  'Y','N'
% Perform the a posteriori detection
% of mixed LDs ?
%
% Default: 'Y;

idpar.mix_detect='y'; 

% idpar.discard_threshold_factor: integer cluster and /or mode datasets
% with less than idpar.discard_threshold_factor * (dimension of PVs) points
% will be discarded. Usually, it is wise to set this parameter equal to 2
% or 3 (as a rule of thumb, if data are noisy, one needs at least 2 or 3
% datapoints for each parameter that must be estimated).
% The minimal value must be 1
%
% Default: 'Y;

idpar.discard_threshold_factor=2; 

% idpar.YNquestions:  'Y','N'
% Ask the user if the regression algorithm must go on after clustering.
% This is useful because if many clusters have been found and one uses
% MRLP, maybe it will take a *very* long time to get the MRLP solution. In
% this case, one would like to be given an option if to continue or not.

idpar.YNquestions='n'; 


                                  
% NEXT, THE STRUCTURE PLOTPAR IS DEFINED.
% CHANGE ONLY THE VALUES BUT *NOT* THE NAME OF THE FIELDS.



% plotpar.marker_modes:  colored markers for plotting classified xi-points
% or mode datasets.
% If the number of markers is less than the number of modes, the markers
% will be replicated (and a warning will be displayed).
% NB: don't use '+' or squares: these markers are reserved for cluster
% centers and for outliers
%

plotpar.marker_modes=['db'; '^b'; 'ob';'*b';'xb';'vb';'>b';'pb';'hb';'<b';...
    'dr'; '^r'; 'or';'*r';'xr';'vr';'>r';'pr';'hr';'<r';...
    'dg'; '^g'; 'og';'*g';'xg';'vg';'>g';'pg';'hg';'<g';];

% plotpar.colors_regions:  RGB triplets for coloring different mode regions
% when the regressor set is 2-dimensional. 
% If the number of triplets is less than the number of modes, the triplets
% will be replicated (and a warning will be displayed).
%

plotpar.color_regions=colorcube(30);

% plotpar.color_surface:  RGB triplet for the color of the
% mode hyperplanes (when the regressor set is 2D) 
%
% Suggested: [0.4 0.4 0.8]. Good for B&W printing.

plotpar.color_surface=[0.4 0.4 0.8];

% plotpar.perm_colors: permute the colors in plotpar.colors_regions this is
% useful for coloring, by trial and error, the regions exactly as in the
% original model (when data are artificially generated).
% To see why this is useful, remember that the identified models may be NOT
% ordered as the original models (for instance the identified mode n. 2 can
% be the what is the original mode n.1).
%
% NB: Obviously, if the original modes are unknown, as in real-world 
% identification, the vector perm_colors does not play any role and
% one can choose perm_colors=1:length(plotpar.colors_regions), i.e.
% no permutation.
%
% Note that if plotpar.perm_colors is not defined by the user, it is
% automatically set to 1:length(plotpar.color_regions)in hit_plot_models

plotpar.perm_colors=1:length(plotpar.color_regions);

%  plotpar.x_grid, plotpar.y_grid: grids for plotting 2D PWA functions. 
%  Empty matrices mean that grids of 30 points will be automatically
%  generated by hit_plot_models and plotpar.x_grid, plotpar.y_grid will be
%  updated accordingly

plotpar.x_grid=[];
plotpar.y_grid=[];

% plotpar.datapoints.yn: 'Y' or 'N'
% Decide if hit_plot_idmodes must also plot the classified datapoints or
% just the modes

plotpar.datapoints_yn='Y';

% Next, set the plots to be displayed during the running of hit_regression

%  plotpar.plot_xi_points.yn: 'Y' or 'N'
%  plotpar.plot_xi_points.fig: figure number
%  These two fields control if xi-points must be displayed before
%  clustering

plotpar.plot_xi_points_yn='N';
plotpar.plot_xi_points_fig=10;

%  plotpar.plot_class_xi_points.yn: 'Y' or 'N'
%  plotpar.plot_class_xi_points.fig: figure number
%  These two fields control if classified xi-points must be displayed after
%  clustering

plotpar.plot_class_xi_points_yn='Y';
plotpar.plot_class_xi_points_fig=11;

printcopyright(hit_ver);
%==========================================================================
function printcopyright(hit_ver)
% prints the copyright notice

fprintf('\nHybrid Identification Toolbox (HIT) %s initialized...\n', hit_ver);
fprintf('Copyright (C) 2005 by G. Ferrari-Trecate');
fprintf('\nEmail: giancarlo.ferrari@unipv.it\n\n');

