function output = callqp_cplex(interfacedata)

% Author Johan Löfberg 
% $Id: callqp_cplex.m,v 1.2 2005/05/07 13:53:20 joloef Exp $

% Retrieve needed data
options = interfacedata.options;
F_struc = interfacedata.F_struc;
c       = interfacedata.c;
Q       = interfacedata.Q;
K       = interfacedata.K;
x0      = interfacedata.x0;
integer_variables = interfacedata.integer_variables;
ub      = interfacedata.ub;
lb      = interfacedata.lb;

showprogress('Calling CPLEXINT',options.showprogress);

% Notation used
C = full(c);                    % Must be full
if ~isempty(F_struc)
    A = (-F_struc(:,2:end));    % Must be full
    B = full(F_struc(:,1));         % Must be full
else
    A = ones(1,length(c)); %FAILS anyway due to bug in CPLEX...
    B = 1e9;
end
CTYPE = repmat('L',K.l+K.f,1);  % Standard variables
CTYPE(1:K.f) = 'E';             % Equality constrained variables

if options.savedebug
    save cplexdebug
end
% Call mex-interface
solvertime = clock; 
try
    [x,FMIN,STATUS,SLACK,DUAL] = qp_cplex(1,2*(Q),C,A,B,CTYPE,lb,ub,[],options.verbose);
    D_struc = [];
catch
    D_struc = [];
    x = nan*ones(length(c),1);
    STATUS = 32201;  
    FMIN = [];
    SLACK = [];
    DUAL = [];
end
solvertime = etime(clock,solvertime);

% Duals available currently
D_struc = [];

% Check, currently not exhaustive...
switch STATUS
case {1,101,102}
    problem = 0;
case {0,2,103}
    problem = 1;
case {3,4}
    problem = 1;
case 119
    problem = 12;
case 32201
    problem = -5;        
otherwise
    problem = -1;
end
infostr = yalmiperror(problem,'CPLEXINT');	

% Save all data sent to solver?
if options.savesolverinput
    solverinput.SENSE = 1;
    solverinput.Q = Q;    
    solverinput.A = A;
    solverinput.C = C;
    solverinput.B = B;
    solverinput.CTYPE = CTYPE;
    solverinput.LB = [];
    solverinput.UB = [];
    solverinput.X0 = [];
    solverinput.SLACK = [];
else
    solverinput = [];
end

% Save all data from the solver?
if options.savesolveroutput
    solveroutput.XMIN = x;
    solveroutput.FMIN = FMIN;
    solveroutput.STATUS = STATUS;
    solveroutput.SLACK = SLACK;
else
    solveroutput = [];
end

% Standard interface 
output.Primal      = x(:);
output.Dual        = D_struc;
output.Slack       = [];
output.problem     = problem;
output.infostr     = infostr;
output.solverinput = solverinput;
output.solveroutput= solveroutput;
output.solvertime  = solvertime;