function [Fhull,t] = hull(varargin)
% HULL  Construct a model of the convex hull
%
% H = hull(F1,F2,...)
%
% OUTPUT
%   H   : SET object describing the convex hull of the input constraints
%
% INPUT
%   Fi  : SET objects with constraints
%
% Note that the convex representation of the convex hull requires a lifting
% (introduction of auxially variables). Hence, if you have many set of
% constraints, your problem rapidly grows large.

% $Id: hull.m,v 1.7 2006/06/21 13:31:03 joloef Exp $   


% Pre-process to convert convex quadratic constraints to socp constraints.
% This makes the perspective code easier
for i = 1:nargin
    varargin{i} = convertquadratics(varargin{i});
end

variables = [];
for i = 1:nargin
    if ~(isa(varargin{i},'lmi') | isa(varargin{i},'socc'))
        error('Hull can only be applied to linear constraints');
    elseif ~(islinear(varargin{i}))
        error('Hull can only be applied to linear and convex quadratic constraints');
    end
    variables = unique([variables depends(varargin{i})]);
end

if nargin == 1
    Fhull = varargin{1};
    return
end

y = sdpvar(repmat(length(variables),1,nargin),repmat(1,1,nargin));
t = sdpvar(nargin,1);

Fhull = set([]);
for i = 1:nargin
    Fi = varargin{i};
    tvariable = getvariables(t(i));
    for j = 1:length(Fi.clauses)
        local_variables = getvariables(Fi);
        Xi = Fi.clauses{j}.data;
        local_variables = getvariables(Xi);
        local_index = find(ismember(variables,local_variables));
        new_variables = getvariables(y{i}(local_index));
        Fi.clauses{j}.data = brutepersp(Fi.clauses{j}.data,tvariable,new_variables);       
    end
    Fhull = Fhull + Fi;
end
Fhull = Fhull + set(sum([y{:}],2) == recover(variables));
Fhull = Fhull + set(sum(t)==1) + set(t>0);