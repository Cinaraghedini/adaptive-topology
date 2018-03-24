function [sep_hyp,spec]=hit_prec_PSVC(X_or,F,s)
%HIT_PREC_PSVC Compute mode regions through Proximal Support Vector Classification (PSVC) 
%
% -------------------------------------------------------------------------
% DESCRIPTION
% -------------------------------------------------------------------------
% [sep_hyp,spec]=hit_prec_PSVC(X_or,F,s)
%
% -------------------------------------------------------------------------
% INPUT
% -------------------------------------------------------------------------
% X_or: matrix containing regressors. Each row is a datapoint.
%
% F: structure containing information about the mode datasets (the
% classified datapoints that are also inliers, i.e. not discarded during
% regression). Type 'help hit_regression' for a description of its fields.
%
% s: number of regions to be found.
%
% -------------------------------------------------------------------------
% OUTPUT                                                                                                   
% -------------------------------------------------------------------------
% sep_hyp{i,j} stores the hyperplane separating the i-th and j-th mode
% datasets according to the formula sep_hyp{i,j}*[x  1]=0. By convention,
% sep_hyp{i,i}=[0 .. 0] thus defining a fictious hyperplane. 
%
% spec: structure containing information about the results
%   spec.correctness(i,j), i>j: correctness in separating regressors of the
%   i-th mode from regressors of the j-th mode. It is the ratio between the
%   number of regressors correctly classified divided the total number of
%   regressors of the i-th and j-th mode.

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

global idpar;
if ~isstruct(idpar),
    hit_error;
end

dim_domain=size(X_or,2);

% ClustBin_i and ClustBin_j are cell arrays used in the subsequent loop 'for'
% 		ClustBin_i{1} contains the x-points belonging to the i-th cluster
% 		ClustBin_i{2} contains the corresponding outputs all equal to 1
%
% 		ClustBin_j{1} contains the x-points belonging to the j-th cluster
% 		ClustBin_j{2} contains the corresponding outputs all equal to -1

ClustBin_i=cell(1,2);
ClustBin_j=cell(1,2);

% spec.correctness stores the training set correctness in separating two
% regions (i.e. the percentage of points separated without errors).
% spec.correctness(i,j) is the correctness of the hyperplane
% separating mode i from mode j.
% Large values are usually caused by many misclassified datapoints.

spec.correctness=zeros(s,s);
for i=1:s

    ClustBin_i{1}=X_or(F.pos{i},:);
    ni=size(F.y{i},1);
    ClustBin_i{2}=ones(ni,1);  %set as 1 the outputs of the i-th cluster

    for j=i+1:s
        %fprintf('\n Computing the hyperplane separating region %1d from region %1d\n',i,j);
        %extract the original x-points belonging to the j-th cluster
        ClustBin_j{1}=X_or(F.pos{j},:);
        nj=size(F.y{j},1);
        ClustBin_j{2}=-ones(nj,1);  %set as -1 the outputs of the j-th cluster
        % in order to use SVM stores all the x-points and the 1/-1 outputs
        % in a single vector
        ClustBin_totx=[ClustBin_i{1};ClustBin_j{1}];
        ClustBin_toty=[ClustBin_i{2};ClustBin_j{2}];

        % Call PSVM: the weight nu must be VERY BIG for having sensible results !
        nu=1e+6;

        [w,gamma, trainCorr, testCorr, cpu_time, nu]=psvc(ClustBin_totx,ClustBin_toty,1,nu,0,1);
        spec.correctness(i,j)=trainCorr;

        % compute the separating hyperplane wij*x+b0=0


        % from psvc.m we have that
        % wij * x-b0 is > 0 on the points in the i-th class.
        % Therefore signs must be inverted to have
        % -wij x+b0 <0 on the i-the class
        %
        %

        wij=-w;
        b0=+gamma;

        %append the displacement to wij

        wij=[wij',b0];


        %store wij in the matrix sep_hyp

        sep_hyp{i,j}=wij;

        % the Hyperplane separating the j-th region from the i-th region has
        % the opposit sign to preserve  wji x +b <0 on the j-th region

        sep_hyp{j,i}=-wij;
    end
end



% add a useless Hyperplane [0 0 .. 0] * x -1 <0
% to fill the (i,i) elements of sep_hyp
% This avoids to complicate the use of indices in the
% determination of the regions with polyreduce
% (if I have time, I will remove this trick)
%
for i=1:s
    wii=zeros(1,dim_domain+1);
    sep_hyp{i,i}=wii ;
end
