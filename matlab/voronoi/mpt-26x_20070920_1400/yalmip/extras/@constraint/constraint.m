function F = constraint(X,quantifier,Y)
% Internal class for constraint list

% Author Johan Löfberg
% $Id: constraint.m,v 1.6 2007/03/12 12:51:53 joloef Exp $

superiorto('sdpvar');
superiorto('double');

if isa(X,'blkvar')
    X = sdpvar(X);
end
if isa(Y,'blkvar')
    Y = sdpvar(Y);
end

% Try to evaluate
switch quantifier
case {'>','>='}
    Z = X - Y;
case {'<','<=','=='}
    Z = Y - X;
otherwise
    error('Quantifier not supported')
end

F.List={X,quantifier,Y};
F.Evaluated{1} = Z;

if isequal(Z,0)
    warning('Constraint evaluated to trivial true.')
    F = set([]);
    return
end

switch quantifier
case {'>','<'}
    F.strict(1) = 1;
case {'>=','<=','=='}
    F.strict(1) = 0;
otherwise
    error('Quantifier not supported')
end

F = class(F,'constraint');
	