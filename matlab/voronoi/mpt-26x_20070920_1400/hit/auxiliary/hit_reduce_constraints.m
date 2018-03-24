function [estimated_regions,adjacences,adjacences_hyp]=hit_reduce_constraints(sep_hyp, Bigregion,solver)
%HIT_REDUCE_CONSTRAINTS Estimate regions from hyperplanes and a bounding polytope.
%
% -------------------------------------------------------------------------
% DESCRIPTION
% -------------------------------------------------------------------------
% [estimated_regions,adjacences,adjacences_hyp]=...
% hit_reduce_constraints(sep_hyp, Bigregion,solver)
%
% -------------------------------------------------------------------------
% INPUT
% -------------------------------------------------------------------------
% sep_hyp{i,j} stores the hyperplane separating the i-th and j-th mode
% datasets according to the formula sep_hyp{i,j}*[x  1]=0. By convention,
% sep_hyp{i,i}=[0 .. 0] thus defining a fictious hyperplane.
% Bigregion: bounding polytope containing the regions (usually the
% regressor set).
% solver: defines the solver to be used in mpt_solveLPs (called by
% polyreduce_ext.m)
%
% -------------------------------------------------------------------------
% OUTPUT
% -------------------------------------------------------------------------
% estimated_regions(i): polytope defining the i-th region.
% adjacences: each row is a pair (i,j), i<j indicating that the regions i
% and j are adjacent The constraint i<j avoids storing both pairs (i,j) and
% (j,i).
% adjacences_hyp: each row contains the vector [w b] of the hyperplane
% separating region i from region j (the regions are stored in the
% corresponding row of adjacences) according to the formula [w b]*[x,1]'=0.

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

% Extract useful dimensions
% mm: # of models
% ndim: dimension of the regressor space
mm=size(sep_hyp,1);
ndim=length(sep_hyp{1,1})-1;

% adjacences: each row is a pair (i,j), i<j indicating that the regions i and j are adjacent
% The constraint i<j avoids storing both pairs (i,j) and (j,i).
adjacences=[];

% adjacences_hyp: each row contains the vector [w b] of the hyperplane
% separating
% region i from region j (the regions are stored in the corresponding row of adjacances)
% according to the formula [w b]*[x,1]'=0
%
adjacences_hyp=[];

% estimated_regions is a cell array. the i-th region is given by
% 		the inequalities
%
% 		  estimated_regions{i,1}*x< estimated_regions{i,2}

%estimated_regions=cell(mm,2);
for i=1:mm
    Ai=zeros(mm,ndim);
    Bi=zeros(mm,1);
    for j=1:mm
        wij=sep_hyp{i,j};
        Ai(j,:)=wij(1:ndim);
        Bi(j)=-wij(ndim+1);
    end
    % add to each region the constraint of the
    % regressor set XX. Note that the fact that the dimension of Bigregion
    % coincide with the dimension of the regressors, has been checked
    % already at the beginning of hit_regression  
    [BigA,BigB]=double(Bigregion);
    Ai=[Ai;BigA];
    Bi=[Bi;BigB];
    %
    % THE NEXT PART COULD BE IMPROVED
    %
    % Indeed, one would like to remove the next line and compute directly
    % estimated_regions(i)=polytope(At,Bt);
    % so avoiding the use of 'polyreduce_ext'.
    % The problem is that 'polytope' removes redundant constraints but it
    % does not keep track of the constraints that survive (to be stored in left_constr) !!
    % That is why polyreduce_ext is still needed
    %

    [At,Bt,how,left_constr]=polyreduce_ext(Ai,Bi,solver);
    switch(lower(how))
        case{'unbounded', 'infeasible'}
            fprintf('\n =========================================================== ');
            fprintf('\n ERROR IN REDUCING THE CONSTRAINTS OF A REGION POLYTOPE:');
            fprintf('\n mpt_solveLPs reported infeasibility or unboundedness');
            fprintf('\n in solving a constraint reduction problem.');
            fprintf('\n hit_regression.m ends here. ');
            fprintf('\n ============================================================= \n');
            error('hit_reduce_constraints:end','hit_reduce_constraints ends here.')
    end
    estimated_regions(i)=polytope(At,Bt);
    % Build adjacences and adjacences_hyp
    %
    % left_constr stores the indices of the constraints that survived.
    % If the lef_constr(j) is <= mm then it is a constraint
    % separating region ii from region left_constr(j) (and not a constraints of idpar.Bigregion).
    % In this case, if (i,j), i<j is not already in adjacences, this pair
    % is added to adjacences and the corresponding hyperplane to
    % adjacences_hyp
    for kk=1:length(left_constr)
        if   left_constr(kk)<= mm % the left constraint do not belong to idpar.Bigregion
            tempp=sort([i,left_constr(kk)]);
            flag=0;
            for ii=1:size(adjacences,1);
                if adjacences(ii,:)==tempp
                    flag=1 ;
                end
            end
            if flag==0
                adjacences=[adjacences;i,left_constr(kk)];
                norm_sep_hyp=sep_hyp{i,left_constr(kk)};
                %stores the hyperplane
                adjacences_hyp=[adjacences_hyp;norm_sep_hyp];
            end

        end
    end % ends for kk ....
end
