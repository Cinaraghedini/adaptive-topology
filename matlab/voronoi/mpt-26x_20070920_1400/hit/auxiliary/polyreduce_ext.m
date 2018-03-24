function [At,Bt,how,left_constr]=polyreduce_ext(A,B,solver)
% POLYREDUCE_EXT Remove redundant constraints in polytopes.
%
%Given a polyhedron Ax<=B, returns an equivalent polyhedron At x<=Bt by
%eliminating redundant constraints
%
%(C) 1999 by A. Bemporad, F. Torrisi, Zurich, 18/2/1999
%
% G. Ferrari-Trecate added among the output the list of the constraints 
% that are not removed (22/11/2000) and made the function compatible with 
% mpt_solveLPs for solving linear programs.


q=size(B,1);
n=size(A,2);
left_constr=[]; %This will store the row- indices of the deleted constraints

if q~=size(A,1),
    error('The A and B matrices must have the same number of rows.')
end

k=1;

%disp('Don''t worry about warning messages')

At=A;
Bt=B;
j=1;
h=0;

% fprintf('Reducing')
while j<=length(Bt)&(length(Bt)>1),
    h=h+1;
    f=At(j,:);
    ii=[1:j-1,j+1:length(Bt)];
    %aux=evalc('[xopt,dummy,how]=lp(-f'',At(ii,:),Bt(ii));');
    if max(abs(At(j,:))) > 1E-5*Bt(j,:)
        %  NOT a row of all zeros
        auxAt=At(ii,:);
        [aux1,aux2]=size(auxAt);
        
        if nargin==3  % the solver string has been specified: use panmip with the
            % appropriate solver
            
            Options.solver               = solver;
            Options.cplex_file.mpsfile   = 'prob.mps';
            Options.cplex_file.runfile   = 'runfile';
            Options.cplex_file.deletefile= 0;
            
            %lb = -1e8*ones(nvar,1);
            %ub= 1e8*ones(nvar,1);
            % CALL THE SOLVER
            
            templb=-1e6*ones(aux1+aux2,1);
            tempub=[1e6*ones(aux2,1);Bt(ii)];
            lb=templb(1:aux2);
            ub=tempub(1:aux2);
            tildalb=templb(aux2+1:length(templb));
            tildaub=tempub(aux2+1:length(tempub));
            
            type=[];
            for ii=1:aux1+aux2, type=[type;'L'];end
            type=[type;type];
            
            vartype=[];
            for ii=1:aux2, vartype=[vartype;'C'];end
            
            [xopt,fopt,lambda,flag,how]=mpt_solveLPs(-f,[eye(aux2);-eye(aux2);auxAt;-auxAt],[ub;-lb;tildaub;-tildalb],[],[],[],solver);
            switch(lower(how))
    case{'unbounded','infeasible'}
        fprintf('\n =========================================================== ');
        fprintf('\n ERROR IN REDUCING REGION CONSTRAINTS:');
        fprintf('\n mpt_solveLPs reported infeasibility or unboundedness.');
        fprintf('\n hit_regression.m ends here. ');
        fprintf('\n ============================================================= \n');
        error('polyreduce_ext:end','polyreduce_ext ends here.')
end
        else % if nargin=3 ... in this case
            %         % the inputs argument are two and the NAF LP solver is used as default
            %           aux=evalc(['[xopt,istate,objlp,clamda,ifail]=e04mbf(' ...
            %             '-1e6*ones(aux1+aux2,1),[1e6*ones(aux2,1);Bt(ii)]' ...
            %             ',zeros(aux2,1),-f,auxAt,0,1e5);']);
            %       
            %       switch ifail
            %       case {0,1,3}
            %          how = 'feasible';
            %       case 2
            %          how = 'unbounded';
            %       case {4,5}
            %          % might also be considered as infeasible %
            %          warning('LP is cycling or too few iterations')
            %          how = 'infeasible';
            %       case {6,7,8}
            %          how = 'infeasible';
            %       otherwise
            %          error('other error code in "ifail" from e04mbf')
            %       end
            %       
            % 
            %           
        end% if nargin=3 ...
        
        val=f*xopt-Bt(j);
    else
        val = 0;
        how = 'feasible';
    end
    
    if strcmp(how,'unbounded')|(val>1e-10),      
        % Lascia il vincolo
        left_constr=[left_constr;h];
        j=j+1;
        %elseif strcmp(how,'infeasible')
        %   error('The polyhedron Ax<=B is empty')
    else 
        % disp(sprintf('Constraint #%d DELETED!',j));
        At(j,:)=[];
        Bt(j)=[];
    end
    % fprintf('.');
end
% fprintf('\n');