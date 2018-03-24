function output = kktqp(interfacedata)
%KKTQP Solver for indefinite QP problems using binary optmization

% Author Johan L�fberg
% $Id: kktqp.m,v 1.1 2006/04/04 19:01:38 joloef Exp $

% ********************************
%% INITIALIZE DIAGNOSTICS IN YALMIP
% ********************************
K = interfacedata.K;
b = interfacedata.F_struc(1+K.f:end,1);
A = -interfacedata.F_struc(1+K.f:end,2:end);
c = interfacedata.c;

n = length(c); 

if ~all(isinf(interfacedata.lb))
    if ~isempty(lb)
        A = [A;-eye(n)];
        b = [b;-interfacedata.lb]
    end
end
if ~all(isinf(interfacedata.ub))
    if ~isempty(ub)
        A = [A;eye(n)];
        b = [b;interfacedata.ub]        
    end
end

% Lazy...
P = polytope(A,b);
[B,L,U] = bounding_box(P);
[A,b] = double(P);

% Formulation here assumes maximization...
Q = -2*interfacedata.Q;
c = -interfacedata.c;

x = sdpvar(n,1); 
bounds(x,L,U);

y = sdpvar(length(b),1);  % Duals 
dy = binvar(length(b),1); % indicater dual ==0 
ds = binvar(length(b),1); % indicater slack==0 
s = b-A*x; % slack 

% Derive bounds on primal slack
[M,m] = derivebounds(s);

% Let us try to derive bounds on the dauls 
% variables. Hmm...
% Ugly, but let's just maximize/minimize linear
% relaxations
F = set(A'*y == Q*x + c) + set(s>0) + set(y>0);%KKT 
F = F + set(s < ds.*M);   % Big M, we know upper bound on s 
F = F + set(dy+ds <= 1);  % Complementary slackness 
F = F + set(0 <= sum(dy) <= n);

% Find dis-joint constraints (silly way...)
for i = 1:length(b)
    j = findrows(abs(A),abs(A(i,:)));
    if length(j)>1
        S{i} = setdiff(j,i);
    else
        S{i} = [];
    end
end

[a1,a2,a3,model] = export(F,-y(i),sdpsettings('relax',1));
solvertime = clock;
for i = 1:length(b)
    if isempty(S{i})
        model_ = model;
        model_.c  = model_.c*0;
        model_.c(n+i)  = -1;
        sol = feval(model.solver.call,model_);
    else
        k = length(S{i});
        Aeq = sparse(1:k,n+S{i},1,k,n+3*length(b));
        beq = zeros(k,1);
        model_ = model;
        model_.F_struc = [beq Aeq;model.F_struc];
        model_.K.f = model.K.f+k;
        model_.c  = model_.c*0;
        model_.c(n+i)  = -1;
        sol = feval(model.solver.call,model_);
    end
    if sol.problem == 0
        My(i,1) =  sol.Primal(n+i);
    else
        My(i,1) =  1e4;
    end
end

F = F + set(y <= dy.*My);

obj = -0.5*(c'*x+b'*y); % ==cost in optimal points 
sol = solvesdp(F,obj);


% **********************************
%% CREATE SOLUTION
% **********************************
output.problem = sol.problem;
output.Primal      = double(x);
output.Dual        = [];
output.Slack       = [];
output.infostr      = yalmiperror(output.problem,'KKTQP');
output.solverinput  = 0;
if interfacedata.options.savesolveroutput
    output.solveroutput.solved_nodes = solved_nodes;
    output.solveroutput.lower = lower;
    output.solveroutput.upper = upper;    
else
    output.solveroutput =[];
end
output.solvertime = etime(clock,solvertime);
