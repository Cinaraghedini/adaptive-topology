function p = polytope(X)
% polytope  Converts set object to polytope object        
%
% P     = polytope(F)
% [P,x] = polytope(F)
%
% P : polytope object (Requires the Multi-parametric Toolbox)
% x : sdpvar object defining the variables in the polytope P.H*x<P.K
% F : set-object with linear inequalities

% Author Johan Löfberg
% $Id: polytope.m,v 1.3 2005/02/04 10:10:27 johanl Exp $


if all(is(X,'element-wise'))% & all(is(X,'linear'))
    f = [];
    for i = 1:length(X)
        if  X.clauses{i}.type==2
            fi =  X.clauses{i}.data;
            f = [f;fi(:)];
        end
    end
    
    B = full(getbase(f));
    p = polytope(-B(:,2:end),B(:,1));
else
    error('polytope can only be applied to SET objects with linear inequalities.')
end