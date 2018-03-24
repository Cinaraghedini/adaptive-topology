function F = eq(X,Y)
% Internal class for constraint lists

% Author Johan Löfberg
% $Id: eq.m,v 1.2 2005/05/01 11:52:01 joloef Exp $

superiorto('sdpvar');
superiorto('double');

if isa(X,'sdpvar') 
    if is(X,'binary')
        F=iff(X,Y);
        return;
    end
end
if isa(Y,'sdpvar') 
    if is(Y,'binary')
        F=iff(Y,X);
        return;
    end
end

error('Equalities can not be used in double-sided constraints')
