function output = callfmincon(model)

% Author Johan Löfberg
% $Id: callfmincon.m,v 1.55 2007/08/21 20:25:15 joloef Exp $

% Retrieve needed data
options = model.options;
F_struc = model.F_struc;
c       = model.c;
K       = model.K;
x0      = model.x0;
Q       = model.Q;
lb      = model.lb;
ub      = model.ub;
monomtable = model.monomtable;

if isempty(model.evaluation_scheme)
    model = build_recursive_scheme(model);
end

model = compress_evaluation_scheme(model);

switch options.verbose
    case 0
        options.fmincon.Display = 'off';
    case 1
        options.fmincon.Display = 'final';
    otherwise
        options.fmincon.Display = 'iter';
end

% Do some pre-calc to be used in calls from fmincon
nonlinearindicies = union(find(model.variabletype~=0),model.evalVariables);
linearindicies    = setdiff(find(model.variabletype==0),nonlinearindicies);
model.nonlinearindicies = nonlinearindicies;
model.linearindicies    = linearindicies;

any_constraints = (K.f+K.l)>0;

model.Anonlinineq = [];
model.bnonlinineq = [];
model.Anonlineq = [];
model.bnonlineq = [];

% Extract linear and nonlinear equality constraints
if K.f>0
    Aeq = -model.F_struc(1:1:K.f,2:end);
    beq = model.F_struc(1:1:model.K.f,1);

    nonlinear_equalities_indicies = find(any(Aeq(:,nonlinearindicies),2));
    model.Anonlineq = Aeq(nonlinear_equalities_indicies,:);
    model.bnonlineq = beq(nonlinear_equalities_indicies);

    Aeq(nonlinear_equalities_indicies,:) = [];
    beq(nonlinear_equalities_indicies,:) = [];
    Aeq(:,nonlinearindicies) = [];
    model.F_struc(1:model.K.f,:) = [];
    model.K.f = 0;
else
    Aeq = [];
    beq = [];
end

% Find nonlinear eualities implied by lower and upper bounds
if ~isempty(ub) & ~isempty(lb)
    nonlinearequality = find(lb(nonlinearindicies) == ub(nonlinearindicies));
    if ~isempty(nonlinearequality)
        for i = 1:length(nonlinearequality)
            model.Anonlineq = [model.Anonlineq;eyev(length(c),nonlinearindicies(nonlinearequality(i)))'];
            model.bnonlineq = [model.bnonlineq;lb(nonlinearindicies(nonlinearequality(i)))];
        end
    end
end

% Extract linear and nonlinear inequality constraints
if model.K.l>0
    A = -model.F_struc(1:model.K.l,2:end);
    b = model.F_struc(1:model.K.l,1);

    nonlinear_inequalities_indicies = find(any(A(:,nonlinearindicies),2));

    model.Anonlinineq = A(nonlinear_inequalities_indicies,:);
    model.bnonlinineq = b(nonlinear_inequalities_indicies);

    A(nonlinear_inequalities_indicies,:) = [];
    b(nonlinear_inequalities_indicies,:) = [];
    A(:,nonlinearindicies) = [];

    model.F_struc(1:model.K.l,:) = [];
    model.K.l = 0;
else
    A = [];
    b = [];
end

if isfield(options.fmincon,'LargeScale')
    if isequal(options.fmincon.LargeScale,'off')
        A = full(A);
        b = full(b);
        Aeq = full(Aeq);
        beq = full(beq);
    end
end

% This helps with robustness in bnb in some cases
x0candidate = zeros(length(c),1);
if ~isempty(lb) & ~isempty(ub)
    bounded = find(~isinf(lb) & ~isinf(ub));
    x0candidate(bounded) = (lb(bounded) + ub(bounded))/2;
    bounded_below = find(~isinf(lb) & isinf(ub));
    x0candidate(bounded_below) = lb(bounded_below) + 0.5;
    bounded_above = find(~isinf(lb) & isinf(ub));
    x0candidate(bounded_above) = lb(bounded_above) + 0.5;
end

if isempty(x0)
    x0 = x0candidate(linearindicies);
else
    if ~isempty(lb) & ~isempty(ub)
        x0((x0 < lb) | (x0 > ub)) = x0candidate((x0 < lb) | (x0 > ub));
    end
    x0 = x0(linearindicies);
end

if ~isempty(lb)
    lb = lb(linearindicies);
end
if ~isempty(ub)
    ub = ub(linearindicies);
end

if options.savedebug
    ops = options.fmincon;
    save fmincondebug model A b Aeq beq x0 lb ub ops
end

showprogress('Calling FMINCON',options.showprogress);

% Precalc for the callbacks
model = setup_fmincon_params(model);
derivative_available = 1;
for i = 1:length(model.evalMap)
    if isempty(model.evalMap{i}.properties.derivative)
        derivative_available = 0;
        break
    end
end
if derivative_available%(model.SimpleQuadraticObjective | model.SimpleNonlinearObjective) & isempty(model.evalMap)
    options.fmincon.GradObj = 'on';    
end
if derivative_available%~isempty(model.evalMap)
    options.fmincon.GradConstr = 'on';   
end

solvertime = clock;
[xout,fmin,flag,output,lambda] = fmincon('fmincon_fun',x0,A,b,Aeq,beq,lb,ub,'fmincon_con',options.fmincon,model);
solvertime = etime(clock,solvertime);

if isempty(nonlinearindicies)
    x = xout(:);
else
    x = zeros(length(c),1);
    for i = 1:length(linearindicies)
        x(linearindicies(i)) = xout(i);
    end
    x = x(1:length(c));
end

problem = 0;

% Internal format for duals
D_struc = [];

% Check, currently not exhaustive...
if flag==0
    problem = 3;
else
    if flag>0
        problem = 0;
    else
        if isempty(x)
            x = repmat(nan,length(c),1);
        end
        if c'*x<-1e10 % Likely unbounded
            problem = 2;
        else          % Probably convergence issues
            problem = 5;
        end
    end
end

% Save all data sent to solver?
if options.savesolverinput
    solverinput.A = A;
    solverinput.b = b;
    solverinput.Aeq = Aq;
    solverinput.beq = beq;
    solverinput.options = options.fmincon;
else
    solverinput = [];
end

% Save all data from the solver?
if options.savesolveroutput
    solveroutput.x = x;
    solveroutput.fmin = fmin;
    solveroutput.flag = flag;
    solveroutput.output=output;
    solveroutput.lambda=lambda;
else
    solveroutput = [];
end

% Standard interface
output = createoutput(x,D_struc,[],problem,'FMINCON',solverinput,solveroutput,solvertime);

function model = compress_evaluation_scheme(model);
scalars = {'exp','log','sin','cos','log2','log10','inverse_internal2'};
for i = 1:length(model.evaluation_scheme)
    if strcmp(model.evaluation_scheme{i}.group,'eval')
        clear fun
        for k = 1:length(scalars)
            for j = 1:length(model.evaluation_scheme{i}.variables)
                fun{k}(j) = strcmp(model.evalMap{model.evaluation_scheme{i}.variables(j)}.fcn,scalars{k});
            end
        end
        for k = 1:length(scalars)
            fun_i = find(fun{k});
            if length(fun_i) > 1
                all_outputs = [];
                all_inputs  = [];
                for j = fun_i
                    all_outputs = [all_outputs model.evalMap{model.evaluation_scheme{i}.variables(j)}.computes];
                    all_inputs  = [all_inputs model.evalMap{model.evaluation_scheme{i}.variables(j)}.variableIndex];
                end
                model.evalMap{model.evaluation_scheme{i}.variables(fun_i(1))}.computes = all_outputs;
                model.evalMap{model.evaluation_scheme{i}.variables(fun_i(1))}.variableIndex = all_inputs;
                model.evaluation_scheme{i}.variables(fun_i(2:end)) = nan;
            end
        end
        model.evaluation_scheme{i}.variables(isnan(model.evaluation_scheme{i}.variables)) = [];
    end
end
