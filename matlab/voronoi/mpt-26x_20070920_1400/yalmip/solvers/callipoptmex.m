function output = callipoptmex(interfacedata)

% Author Johan Löfberg
% $Id: callipoptmex.m,v 1.4 2006/03/31 07:25:25 joloef Exp $

% Retrieve needed data
options = interfacedata.options;
F_struc = interfacedata.F_struc;
c       = interfacedata.c;
K       = interfacedata.K;
x0      = interfacedata.x0;
Q       = interfacedata.Q;
lb      = interfacedata.lb;
ub      = interfacedata.ub;
monomtable = interfacedata.monomtable;

% Do some pre-calc to be used in calls from ipopt
temp = sum(interfacedata.monomtable,2)>1;
temp = (sum(interfacedata.monomtable,2)==1 & sum(interfacedata.monomtable~=0,2)==1);
nonlinearindicies = find(full(~temp));
linearindicies = find(full(temp));

% Extract gradient info
%Gradients = F_struc(K.l+K.f+1:1:end,:);
%F_struc(K.l+K.f+1:end,:) = [];

% Bounds explictly supplied?
if isempty(lb) | all(isinf(lb))
    %No, extract from model
    if K.l > 0
        b = F_struc(K.f+1:K.l+K.f,1);
        A = -F_struc(K.f+1:K.l+K.f,2:end);
        [lower,upper,bound_rows] = find_variable_bounds(A,b,[],[]);
        lb(linearindicies) = lower(linearindicies);
        ub(linearindicies) = upper(linearindicies);
        bound_rows = bound_rows(find(~any(A(bound_rows,nonlinearindicies),2)));
        F_struc(K.f + bound_rows,:) = [];
        K.l = K.l - length(bound_rows);
    end
end

% Introduce new variables to handle additional inequalties
originalNvariables = length(c);
originalLinear = linearindicies;
nnew = 0;
if K.l
    F_struc = [F_struc [zeros(K.f,K.l);-eye(K.l)]];
    lb(originalNvariables+1:originalNvariables+K.l) = 0;
    ub(originalNvariables+1:originalNvariables+K.l) = inf;
    c(originalNvariables+1:originalNvariables+K.l) = 0;
    Q = blkdiag(Q,spalloc(K.l,K.l,0));
    monomtable = blkdiag(monomtable,eye(K.l));
    linearindicies = [linearindicies(:)' (originalNvariables+1):(originalNvariables+K.l)]';
    nnew = K.l;
    K.f = K.f + K.l;
    K.l = 0;
end
newinterfacedata = interfacedata;
newinterfacedata.monomtable = monomtable;
newinterfacedata.c = c;
newinterfacedata.Q = Q;
newinterfacedata.linearindicies = linearindicies;
newinterfacedata.nonlinearindicies = nonlinearindicies;
newinterfacedata.c = c;
newinterfacedata.K = K;
newinterfacedata.F_struc = F_struc;



% Ok, now we have a model in IPOPTs prefered format

if isempty(x0)
    x0 = [];
%    x0 = zeros(length(linearindicies),1);
else
    x0 = [];
%    x0 = x0(linearindicies);
end

if ~isempty(lb)
    lb = lb(linearindicies);
end
if ~isempty(ub)
    ub = ub(linearindicies);
end

if options.savedebug
    ops = options.ipopt;
    save ipoptdebug interfacedata A b Aeq beq x0 lb ub ops
end

functiondata.nonlinearindicies = nonlinearindicies;
functiondata.linearindicies = linearindicies;
functiondata.monomtable = monomtable;
functionGradientdata.nonlinearindicies = nonlinearindicies;
functionGradientdata.linearindicies = linearindicies;
functionGradientdata.monomtable = monomtable;
functionHessiandata.nonlinearindicies = nonlinearindicies;
functionHessiandata.linearindicies = linearindicies;
functionHessiandata.monomtable = monomtable;

constraintdata.nonlinearindicies = nonlinearindicies;
constraintdata.linearindicies = linearindicies;
constraintdata.monomtable = monomtable;
constraintGradientdata.nonlinearindicies = nonlinearindicies;
constraintGradientdata.linearindicies = linearindicies;
constraintGradientdata.monomtable = monomtable;
constraintHessiandata.nonlinearindicies = nonlinearindicies;
constraintHessiandata.linearindicies = linearindicies;
constraintHessiandata.monomtable = monomtable;

functiondata.F_struc = [interfacedata.f c(:)'];
constraintdata.F_struc = F_struc(1:K.l+K.f,:);

start = 1;
start = start+originalNvariables;
datasaver(functiondata,[],[],constraintdata,[],[])

i_lb = find(~isinf(lb));
val_lb = lb(i_lb);
i_ub = find(~isinf(ub));
val_ub = ub(i_ub);

NoEquCons = K.f;
StartPoint = zeros(length(lb),1);

solvertime = clock;
switch options.verbose
    case 1
%        [xout,ObjFinal,mlb,mub,lambda,info,XInfo,ObjInfo] = Ipopt('ipopt_fun','ipopt_con','fd','fd','fd', NoEquCons, StartPoint, i_lb, val_lb,  i_ub,val_ub,setup_fmincon_params(newinterfacedata));%,'print',-1,'IPRINT',-1);
        [xout,ObjFinal,mlb,mub,lambda,info,XInfo,ObjInfo] = Ipopt('ipopt_fun','ipopt_con','fmincon_fungrad','fd','fmincon_funhessian', NoEquCons, StartPoint, i_lb, val_lb,  i_ub,val_ub,setup_fmincon_params(newinterfacedata));%,'print',1,'IPRINT',0);
    otherwise
        [xout,ObjFinal,mlb,mub,lambda,info,XInfo,ObjInfo] = Ipopt('ipopt_fun','ipopt_con','fmincon_fungrad','fd','fmincon_funhessian', NoEquCons, StartPoint, i_lb, val_lb,  i_ub,val_ub,setup_fmincon_params(newinterfacedata));%,'print',1,'IPRINT',0);
%         [xout,ObjFinal,mlb,mub,lambda,info,XInfo,ObjInfo] =
%         Ipopt('ipopt_fun','ipopt_con','fd','fd','fd', NoEquCons, StartPoint, i_lb, val_lb,  i_ub,val_ub,setup_fmincon_params(newinterfacedata));
end
        
solvertime = etime(clock,solvertime);

% get rid of slaxks
xout = xout(1:length(originalLinear));
linearindicies = originalLinear;

% Put the linear variables in place
if isempty(nonlinearindicies)
    x = xout(:);
else
    x = zeros(length(interfacedata.c),1);
    for i = 1:length(linearindicies)
        x(linearindicies(i)) = xout(i);
    end
end

problem = 0;

% Internal format for duals
D_struc = [];

infostr = yalmiperror(problem,'IPOPT');

% Save all data sent to solver?
if options.savesolverinput

else
    solverinput = [];
end

% Save all data from the solver?
if options.savesolveroutput
    solveroutput.xout = x;
    solveroutput.ObjFinal = ObjFinal;
    solveroutput.mlb = mlb;
    solveroutput.mub=mub;
    solveroutput.lmbda=lambda;
    solveroutput.info=info;
    solveroutput.XInfo=ObjInfo;
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

% Some pre-calc that is sent to function evaluations etc
function constant_data = setup_fmincon_params(interfacedata)

monomtable = interfacedata.monomtable;
nonlinearindicies = interfacedata.nonlinearindicies;
linearindicies = interfacedata.linearindicies;

% actuallyused = any(interfacedata.F_struc(:,2:end),1);
% if ~isempty(interfacedata.Anonlinineq)
%     actuallyused =  actuallyused | any(interfacedata.Anonlinineq,1);
% end
% if ~isempty(interfacedata.Anonlineq)
%     actuallyused =  actuallyused | any(interfacedata.Anonlineq,1);
% end
% nonlinearindicies = intersect(nonlinearindicies,find(actuallyused));

constant_data.interfacedata = interfacedata;
constant_data.monomtable = monomtable;
constant_data.nonlinearindicies = nonlinearindicies;
constant_data.linearindicies = linearindicies;


