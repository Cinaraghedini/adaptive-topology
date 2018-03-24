function output = callcplexint(interfacedata)

% Author Johan Löfberg 
% $Id: callcplexint.m,v 1.19 2006/07/13 14:48:55 joloef Exp $

% Retrieve needed data
options = interfacedata.options;
F_struc = interfacedata.F_struc;
c       = interfacedata.c;
Q       = interfacedata.Q;
K       = interfacedata.K;
x0      = interfacedata.x0;
integer_variables = interfacedata.integer_variables;
binary_variables = interfacedata.binary_variables;
UB      = interfacedata.ub;
LB      = interfacedata.lb;

showprogress('Calling CPLEXINT',options.showprogress);

SENSE = 1;     
C = full(c);   
if K.l+K.f+K.q == 0
    A = zeros(1,length(c));A(1)=1;
    B = 1e6;
 %   A = [];
 %   b = [];
else
    A =-(F_struc(1:K.f+K.l,2:end)); 
    B = full(F_struc(1:K.f+K.l,1));            
end

INDEQ = [];
if K.f>0
    INDEQ(1:K.f) = 1:K.f;
end

if K.q(1)>0
    top = K.f+K.l + 1;
    for i = 1:length(K.q)
        % [cx+d;Ax+b]   |Ax+b|<cx+d, originally a QCQP
        m = K.q(i);
        ci = F_struc(top,2:end)';
        di = F_struc(top,1);
        Ai = F_struc(top+1:top+m-1,2:end);
        bi = F_struc(top+1:top+m-1,1);        
        QC(i).Q = full(Ai'*Ai - ci*ci');
        QC(i).r = full(di'*di - bi'*bi);   
        QC(i).L = full(2*bi'*Ai - 2*di*ci');
        top = top+m;
    end
else
    QC = [];
end

VARTYPE = repmat('C',size(A,2),1);
VARTYPE(integer_variables)='I'; % Integer variables
VARTYPE(binary_variables) ='B'; % Binary variables

if nnz(Q)==0
    H = [];
else
    H = full(2*Q);
end

PARAM = options.cplex.param;
OPTIONS.verbose = options.verbose;
OPTIONS.logfile = options.cplex.logfile;
if ~isempty(x0)
    OPTIONS.x0 = [(1:length(x0))' x0(:)];
end

if options.savedebug
    save cplexintdebug H C A B LB UB QC VARTYPE INDEQ PARAM OPTIONS
end

% Call mex-interface
solvertime = clock; 
[x,FMIN,SOLSTAT,DETAILS] = cplexint(H, C, A, B, INDEQ, QC, LB, UB,VARTYPE,PARAM,OPTIONS);
if interfacedata.getsolvertime solvertime = etime(clock,solvertime);else solvertime = 0;end
problem = 0;
D_struc = -DETAILS.dual;    

switch SOLSTAT
    case {1,101,102}
        problem = 0;
    case {3,22}
        problem = 1;
    case 108
        problem = 3;
    case {2}
        problem = 2;
    case {4,119}
        problem = 12;
    otherwise
        problem = -1;
end
infostr = yalmiperror(problem,'CPLEXINT');	

% Save all data sent to solver?
if options.savesolverinput
    solverinput.H = H;
    solverinput.A = A;
    solverinput.C = C;
    solverinput.INDEQ = INDEQ;
    solverinput.QC = QC;
    solverinput.B = B;
    solverinput.CTYPE = CTYPE;
    solverinput.LB = LB;
    solverinput.UB = UB;
else
    solverinput = [];
end

% Save all data from the solver?
if options.savesolveroutput
    solveroutput.x = x;
    solveroutput.FMIN = FMIN;
    solveroutput.SOLSTAT = SOLSTAT;
    solveroutput.DETAILS=DETAILS;
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