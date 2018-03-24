function [final_models_cont]=hit_compute_cont_modes(F,ndata,s,weig,adjacences,adjacences_hyp)
%HIT_COMPUTE_CONT_MODES Estimate the mode PVs in continuous PWA models 
%
% -------------------------------------------------------------------------
% DESCRIPTION
% -------------------------------------------------------------------------
% [final_models_cont]=...
% hit_compute_cont_modes(F,ndata,s,weig,adjacences,adjacences_hyp)
% -------------------------------------------------------------------------
% INPUT
% -------------------------------------------------------------------------
% F: structure containing information about the mode datasets, see
% hit_regression.m
% ndata: total number of datapoints.
% s: number of modes to find.
% weig(i): weight of the i-th datapoint used in weighted LS for estimating
% the mode PVs.
% adjacences: each row is a pair (i,j), i<j indicating that the regions i
% and j are adjacent The constraint i<j avoids storing both pairs (i,j) and
% (j,i).
% adjacences_hyp: each row contains the vector [w b] of the hyperplane
% separating region i from region j (the regions are stored in the
% corresponding row of adjacences) according to the formula [w b]*[x,1]'=0.
%
% -------------------------------------------------------------------------
% OUTPUT                                                                                                   
% -------------------------------------------------------------------------
% final_models_cont(i,:) estimated PV of the i-th mode. 

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

global idpar ;
if ~isstruct(idpar),
    hit_error;
end

% number of modes
m=size(F.X,1);

if ~isempty(adjacences_hyp) %this is the case when more than one mode has been found
    % dimension of regressors
    ndim=size(adjacences_hyp(1,:),2)-1;
    disp('FIND THE CONTINUOS PWARX MODEL')
    % Create the constraints for the adjacent regions
    [Acons,Bcons,Aeq,Beq]=hit_create_constr(s,adjacences,adjacences_hyp);
    % Construct the hessian
    PHI=zeros(ndata,s*(ndim+1));
    PHIcl=0; %count the added lines to PHI
    WEIGHT=[];
    Y=[];
    for i=1:m
        nC=size(F.X{i},1);
        if nC <= ndim+1
            fprintf('WARNING: the data of mode %g are less than the parameters \n',i);
        end
        Phitemp=[F.X{i},ones(nC,1)];
        PHI(PHIcl+1:PHIcl+nC,(1+(i-1)*(ndim+1)):i*(ndim+1))=Phitemp;
        weights=1./weig(F.pos{i}); %weights equal to the confidence measures
        %Weight_primal(F{i,4})=1./weights;
        W=diag(1./weights);

        %  uncomment the next line to NOT use the weights, i.e. set them equal to 1
        %W=eye(nC);
        WEIGHT=[WEIGHT,zeros(size(WEIGHT,1),size(W,2));zeros(size(W,1),size(WEIGHT,2)),W];
        Y=[Y;F.y{i}];
        PHIcl=PHIcl+nC;
    end

    % Compute the QP to find the coefficients
    HESSIAN=PHI'*WEIGHT*PHI;
    LINEAR=-Y'*WEIGHT*PHI;



    % CALL THE SOLVER

    % Syntax for mpt_solveqp
    % With rescue=1, mpt tries different QP solvers, if one reports
    % infeasibility or unboundedness
    rescue=1;
    [Solution,lambda,how,exitflag,FVAL]=mpt_solveQP(HESSIAN,LINEAR,[],[],Aeq,Beq,[],idpar.QPsolver,[],rescue);
    switch(lower(how))
        case{'unbounded', 'infeasible'}
            fprintf('\n =========================================================== ');
            fprintf('\n ERROR IN COMPUTING THE CONTINUOUS PWA MODEL:');
            fprintf('\n mpt_solveQP reported infeasibility or unboundedness');
            fprintf('\n Try to change the QP solver (modifying idpar.QPsolver)');
            fprintf('\n hit_regression.m ends here. ');
            fprintf('\n ============================================================= \n');
            error('hit_regression:end','hit_regression ends here.')
    end

    % Extract the models from the Solution vector and store them in
    % final_models_cont

    final_models_cont=zeros(s,(ndim+1));
    for k=1:s
        temppp=Solution(1+(k-1)*(ndim+1):k*(ndim+1));
        final_models_cont(k,:)=temppp';
    end


else % only one mode has been found
    nC=size(F.X{1},1);
    ndim=size(mX,2);
    PHI=[F.X{1},ones(nC,1)];
    WEIGHT=diag(1./weig(F.pos{1})); % weights equal to the scalar confidence measures
    Y=[F.y{1}];
    HESSIAN=PHI'*WEIGHT*PHI;
    LINEAR=PHI'*WEIGHT*Y;
    temppp=HESSIAN\LINEAR;
    final_models_cont=temppp';
end




%==========================================================================
function [Acons,Bcons,Aeq,Beq]=hit_create_constr(s,adjacences,adjacences_hyp)

% Gives back the equality constraints in the inequality form Acons*weights<=Bcons
% and in the equality form Aeq*weights=Beq. The latter can be used with quadprog.m
global idpar;
if ~isstruct(idpar),
    hit_error;
end

nw=size(adjacences_hyp,2);
nwbar=size(adjacences_hyp,2)-2; % dimension of the vector wbar

%if nwbar<1
%    warning('Dimension of the domain less than 2 ! Create_constr.m cannot be used');
%end

nhyp=size(adjacences,1); % Extract the # of adjacent hyperplanes

Acons=[];
Bcons=[];
Aeq=[];
Beq=[];

for i=1:nhyp
    sep_hyp=adjacences_hyp(i,:); % row vector
    % normalize with respect the first coefficient: IT IS SUPPOSED TO BE DIFFERENT FROM ZERO
    if sep_hyp(1)==0
        warning('FIRST COEFF OF THE HYPERPLANE ==0:  hit_compute_cont_modes.m (Create_constr) reset the coeff to 1e-8');
        sep_hyp(1)=1e-8;
    end

    sep_hyp=sep_hyp/sep_hyp(1);

    Atemp=zeros(2+2*nwbar,s*nw);
    Btemp=zeros(2+2*nwbar,1);
    Aeqtemp=zeros(1+nwbar,s*nw);
    Beqtemp=zeros(1+nwbar,1);

    temp=[-sep_hyp(nw), zeros(1,nwbar),1];
    tempp= [-(sep_hyp(2:nwbar+1))',eye(nwbar),zeros(nwbar,1)];

    %positioning of the constraints on the k-th hyperplane
    pos=adjacences(i,1);
    Atemp(1,1+(pos-1)*nw:nw*pos)=-temp;
    Atemp(2,1+(pos-1)*nw:nw*pos)=temp;
    Atemp(3:2+nwbar,1+(pos-1)*nw:nw*pos)=-tempp;
    Atemp(3+nwbar:2+2*nwbar,1+(pos-1)*nw:nw*pos)=tempp;

    Aeqtemp(1,1+(pos-1)*nw:nw*pos)=-temp;
    Aeqtemp(2:1+nwbar,1+(pos-1)*nw:nw*pos)=-tempp;

    %positioning of the constraints on the j-th hyperplane
    pos=adjacences(i,2);
    Atemp(1,1+(pos-1)*nw:nw*pos)=temp;
    Atemp(2,1+(pos-1)*nw:nw*pos)=-temp;
    Atemp(3:2+nwbar,1+(pos-1)*nw:nw*pos)=tempp;
    Atemp(3+nwbar:2+2*nwbar,1+(pos-1)*nw:nw*pos)=-tempp;

    Aeqtemp(1,1+(pos-1)*nw:nw*pos)=temp;
    Aeqtemp(2:1+nwbar,1+(pos-1)*nw:nw*pos)=tempp;

    Acons=[Acons;Atemp];
    Bcons=[Bcons;Btemp];

    Aeq=[Aeq;Aeqtemp];
    Beq=[Beq;Beqtemp];

end % ends for i=...
