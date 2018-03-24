function indicies = find(x)
%FIND (overloaded)

% Author Johan Löfberg 
% $Id: find.m,v 1.3 2005/10/18 16:44:27 joloef Exp $   

base = x.basis;
vars = x.lmi_variables;
indicies = find(any(base,2));
