function [sep_hyp,spec]=hit_prec_MRLP(X_or,F,s);
%HIT_PREC_MRLP Compute mode regions through Multicategory Robust Linear Programming (MRLP) 
%
% -------------------------------------------------------------------------
% DESCRIPTION
% -------------------------------------------------------------------------
% [sep_hyp,spec]=hit_prec_MRLP(X_or,F,s);
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


% Build the constraints of the LP
%
%				min f*x
%
%				s.t. A*x<=b
%


A=[];
b=[];
f=[];


dim=size(X_or,2);

m=[];
for i=1:s
    m=[m;size(F.pos{i},2)];
end
dimY=0; %compute the lenght of the vector of the slack vaiables Y
for i=1:s
    dimY=dimY+m(i)*(s-1);
end

for i=1:s
    P=X_or(F.pos{i},:);
    for j=1:s
        if i~=j
            count=0;
            for k=1:i-1
                count=count+(s-1)*m(k);
            end
            if j<i
                count=count+(j-1)*(m(i));
            else
                count=count+(j-2)*(m(i));
            end

            Atemp=zeros(m(i),s*dim+dimY+s);
            btemp=-ones(m(i),1);
            Atemp(:,dim*(i-1)+1:dim*i)=-P;
            Atemp(:,dim*(j-1)+1:dim*j)=P;
            Atemp(:,dim*s+count+1:dim*s+count+m(i))=-eye(m(i));
            Atemp(:,dim*s+dimY+i)=ones(m(i),1);
            Atemp(:,dim*s+dimY+j)=-ones(m(i),1);
            % at 06/12/05, sparse matrices are not handled very well by 
            % the LP solver of NAG. Then, the following line is commented
            %
            % A=sparse([A;Atemp]);
            A=[A;Atemp];
            b=[b;btemp];
        end %ends if i~=j
    end %ends for j
end %ends for

%positivity constraints for yij

Atemp=[zeros(dimY,s*dim),-eye(dimY),zeros(dimY,s)];

btemp=zeros(dimY,1);
A=[A;Atemp];
b=[b;btemp];
if exist('count')
    fprintf('\n Number of constraints in the MRLP problem : %f\n',count);
end
% cost functional
tempf=[];
for i=1:s
    tempf=[tempf,ones(1,m(i)*(s-1))/m(i)];
end

f=[zeros(1,s*dim),tempf,zeros(1,s)];
nvar=size(f,2);
nconstr=size(b,1);
%%New part
% ctype=[];
% for i=1:nconstr
%     ctype=[ctype;'L'];
% end
% 
% vartype=[];
% for i=1:nvar
%     vartype=[vartype;'C'];
% end


% CALL THE LP SOLVER
% fprintf('MRLP: Call the LP solver for finding the regions\n');

[xopt,fopt,lambda,flag,how]=mpt_solveLPs(f,A,b,[],[],[],idpar.LPsolver);

switch(lower(how))
    case{'unbounded','infeasible'}
        fprintf('\n =========================================================== ');
        fprintf('\n ERROR IN COMPUTING ONE HYPERPLANE SEPARATING THE REGIONS:');
        fprintf('\n mpt_solveLPs reported infeasibility or unboundedness');
        fprintf('\n in solving an MRLP problem.');
        fprintf('\n Suggestion: try with a different pattern recognititon');
        fprintf('\n algorithms (e.g. SVC or PSVC).');
        fprintf('\n hit_regression.m ends here. ');
        fprintf('\n ============================================================= \n');
        error('hit_prec_MRLP:end','hit_prec_MRLP ends here.')
end
% end NEWPART

% aux=evalc('[xopt,fval,exitflag]=linprog(f'',A,b);')

% Extract the wj (store them in the rows of W)
% and the gj (store them in the rows of G)

W=[];
G=[];
for i=1:s
    W=[W;(xopt((i-1)*dim+1:(i)*dim))'];
    G=[G;xopt(dim*s+dimY+i)];
end

% stores the hyperplanes
% sep_hyp: cell array. the ij-entry is the hyperplane [w b]*[x,1]'=0
% separating the cluster i from the cluster j
%

% Note that the separating hyperplanes must be estimated on the basis of
% the NON normalized data !!

sep_hyp=cell(s,s);


for i=1:s
    for j=i+1:s
        wij=[W(j,:)-W(i,:),G(i)-G(j)];
        sep_hyp{i,j}=wij;
        sep_hyp{j,i}=-wij;
    end
end


% add an useless Hyperplane [0 0 .. 0] * x -1 <0
% to fill the (i,i) elements of seo_hyp
% This avoids to complicate the use of indices in the
% determination of the regions with polyreduce
% (if I will have time, I will remove this trick)
%


for i=1:s
    wii=zeros(1,dim+1);
    sep_hyp{i,i}=wii ;
end

% spec.correctness stores the training set correctness in separating two
% regions (i.e. the percentage of points separated without errors).
% spec.correctness(i,j) is the correctness of the hyperplane
% separating mode i from mode j.
% Large values are usually caused by many misclassified datapoints.

spec.correctness=zeros(s,s);

for i=1:s
    ClustBin_i{1}=X_or(F.pos{i},:);
    ni=size(F.y{i},1);
    ClustBin_i{2}=ones(ni,1);  
    for j=i+1:s
        ClustBin_j{1}=X_or(F.pos{j},:);
        nj=size(F.y{j},1);
        ClustBin_j{2}=-ones(nj,1);  
        ClustBin_totx=[ClustBin_i{1};ClustBin_j{1}];
        ClustBin_toty=[ClustBin_i{2};ClustBin_j{2}];
        spec.correctness(i,j)=correctness(ClustBin_totx,ClustBin_toty,sep_hyp{i,j}(1:dim)',sep_hyp{i,j}(dim+1));
    end
end
%%%%%%%%%%%%%%% correctness calculation %%%%%%%%%%%%%%%%%%%%
function corr = correctness(AA,dd,w,gamma)

p=-sign(AA*w+gamma);
corr=length(find(p==dd))/size(AA,1)*100;

return
