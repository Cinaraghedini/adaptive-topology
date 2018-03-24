function output = calllpsolve(interfacedata)

% Author Johan Löfberg 
% $Id: calllpsolve.m,v 1.9 2006/05/23 14:15:22 joloef Exp $

% Retrieve needed data
options = interfacedata.options;
F_struc = interfacedata.F_struc;
c       = interfacedata.c;
K       = interfacedata.K;
integer_variables = interfacedata.integer_variables;
binary_variables = interfacedata.binary_variables;
ub      = interfacedata.ub;
lb      = interfacedata.lb;

n = length(c);
% Bounded variables converted to constraints
if ~isempty(ub)
    LB = lb;
    UB = ub;
   % LB(isinf(LB)) = -1e12;
   % UB(isinf(UB)) = 1e12;
else
    if isempty(integer_variables) & isempty(binary_variables)         
        LB = -ones(n,1)*inf;;%-ones(n,1)*0;
        UB = ones(n,1)*inf;;%[];%ones(n,1)*500;
    else
        %LP_SOLVE FAILS IF BOUNDS NOT EXPLICIT
        [LB,UB,used_rows] = findulb(F_struc,K);
        LB(isinf(LB)) = -1e12;
        UB(isinf(UB)) = 1e12;
        F_struc(K.f+used_rows,:)=[];
        K.l = K.l - length(used_rows);
        LB(binary_variables) = max(LB(binary_variables),0);
        UB(binary_variables) = min(UB(binary_variables),0);
    end
end

if options.showprogress;showprogress('Calling LPSOLVE',options.showprogress);end

f = - full(c);            % Must be full
A = - F_struc(:,2:end);
b = full(F_struc(:,1));   % Must be full
e = -ones(size(A,1),1);
e(1:K.f) = 0;

xint = uniquestripped([integer_variables binary_variables]);

% Call mex-interface
solvertime = clock; 
if options.savedebug
    save mxlpsolvedebug f A b e UB LB xint 
end

lp = create_lp_solve_model(A,b,f,xint,LB,UB,e,options);

try
    solvertime = clock; 
    if options.showprogress;showprogress(['Calling ' interfacedata.solver.tag],options.showprogress);end
    result=mxlpsolve('solve', lp);
    solvertime = etime(clock,solvertime);
    if result == 0 | result == 1 | result == 11 | result == 12
        [obj, x, duals] = mxlpsolve('get_solution', lp);
    else
        obj = [];
        x = zeros(length(c),1);
        duals = [];
    end
    mxlpsolve('delete_lp', lp);
catch
    obj = [];
    x = zeros(length(c),1);
    duals = [];
    result = -1;
    mxlpsolve('delete_lp', lp);
end

if options.saveduals & isempty(integer_variables)
    D_struc = duals;
else
    D_struc = [];
end

switch result
    case 0
        problem = 0; % OPTIMAL
    case 2
        problem = 1; % INFEASIBLE
    case 3
        problem = 2; % UNBOUNDED
    case {1,7,12,13}
        problem = 3; % RUN OUT OF TIME OR SIMILIAR
    case 5
        problem = 4;
    case {-2,10,11}
        problem = 11;
    otherwise
        problem = -1;
end
        
infostr = yalmiperror(problem,interfacedata.solver.tag);

% Save all data sent to solver?
if options.savesolverinput
	solverinput.A = A;
	solverinput.f = f;
	solverinput.b = b;	
	solverinput.LB = LB;
	solverinput.UB = UB;    
    solverinput.xint = xint;    
else
	solverinput = [];
end

% Save all data from the solver?
if options.savesolveroutput
	solveroutput.x = x;
    solveroutput.obj = obj;
    solveroutput.duals = duals;  
    solveroutput.result = result;  
else
	solveroutput = [];
end


% Standard interface 
output.Primal      = x;
output.Dual        = D_struc;
output.Slack       = [];
output.problem     = problem;
output.infostr     = infostr;
output.solverinput = solverinput;
output.solveroutput= solveroutput;
output.solvertime  = solvertime;