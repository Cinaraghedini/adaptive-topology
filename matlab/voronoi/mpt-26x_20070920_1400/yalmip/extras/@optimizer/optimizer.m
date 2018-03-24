function sys = optimizer(F,h,ops,x,u)
%OPTIMIZER  Container for optimization problem
%
%   [OPT,PROBLEM] = OPTIMIZER(F,h,options,x,u) exports an object that
%   contains precompiled numerical data to be solved for varying arguments
%   x, returning the optimal value of the expression u.
%
%   OPTIMIZER typically only makes sense if the varying data x enters the
%   optmization problem affinely.
%
%   Example
%    The following problem creates an LP with varying upper and lower
%    bounds on the decision variable.
%
%    The optimizing argument is obtained by indexing (with {}) the optimizer 
%    object with the point of interest. The argument should be a column
%    vector (if the argument has a width larger than 1, YALMIP assumes that
%    the optimal solution should be computed in several points) 
%   
%     A = randn(10,3);
%     b = rand(10,1)*19;
%     c = randn(3,1);
%
%     z = sdpvar(3,1);
%     sdpvar UB LB
%
%     F = set(A*z < b) + set(LB < z < UB);
%     h = c'*z
%     optZ = optimizer(F,h,[],[LB; UB],z);
%     
%     % Solve for LB=1, UB = 3;
%     zopt = optZ{[1; 3]}
%
%     % Solve for (LB,UB) [1;3] and (LB,UB) [2;6]
%     zopt = optZ{[[1; 3], [2;6]]}

if nargin < 5
    error('OPTIMIZER requires 5 inputs');
end

[n,m] = size(x);
x = x(:);

if isempty(ops)
    ops = sdpsettings;
end
ops.avoidequalitybounds=1;

if isa(F,'constraint')
    F = set(F);
end

if any(is(F,'uncertain'))
    [F,h,failure] = robustify(F,h,ops);
    [aux1,aux2,aux3,model] = export(set(x == repmat(pi,n*m,1))+F,h,ops,[],[],0);
else
    [aux1,aux2,aux3,model] = export(set(x == repmat(pi,n*m,1))+F,h,ops,[],[],0);    
end

if norm(model.F_struc(1:n,1)-repmat(pi,length(x),1),inf) > 1e-10
    error('Failed exporting the model (try to specify another solver)')    
end

% Try to set up an optimal way to compute the output
base = getbase(u);
if is(u,'linear') & all(sum(base | base,2) == 1) & all(sum(base,2))==1 & all(base(:,1)==0)
    % This is just a vecotr of variables
    z = [];
    map = [];
    for i = 1:length(u)
        var = getvariables(u(i));
        map = [map;find(var == model.used_variables)];
    end
else
    % Some expression which we will use assign and double to evaluate
    z = recover(depends(u));
    map = [];
    for i = 1:length(z)
        var = getvariables(z(i));
        map = [map;find(var == model.used_variables)];
    end        
end

sys.recover = aux2;
sys.model = model;
sys.dimin = [n m];
sys.dimout = size(u);
sys.map = map;
sys.input.expression = x;
sys.output.expression = u;
sys.output.z = z;
[a,b,c] = find(sys.model.F_struc(1:prod(sys.dimin),2:end));
sys.parameters = b;
used_in = find(any(sys.model.monomtable(:,b),2));
if any(sum(sys.model.monomtable(used_in,:) | sys.model.monomtable(used_in,:),2) > 1)
    sys.nonlinear = 1;
else
    sys.nonlinear = 0;
end
sys = class(sys,'optimizer');

