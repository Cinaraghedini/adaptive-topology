function [properties,F,arguments,fcn]=model(X,method,options,extstruct)
%MODEL  Extracts nonlinear operator models
%
% [properties,F] = model(x)
%
% MODEL returns the constraints needed to model a variable related to an
% extended operator such as min, max, abs, norm, geomean, ...
%
% Examples :
%
% sdpvar x y;
% t = min(x,y);
% [properties,F] = epigraph(t)
% Gives (F = set(t<x) + set(t<y))
%
% sdpvar x y
% t = max(norm([x;y],1+y))
% [properties,F] = epigraph(t)
% Gives (F = set(u<t) + set(1+y<t))
% where u is the variable modelling norm([x;y])

% Author Johan Löfberg
% $Id: model.m,v 1.60 2007/08/09 15:36:37 joloef Exp $


extvar = getvariables(X);
arguments   = cell(1,length(extvar));
properties  = cell(1,length(extvar));

if nargin<2
    method = 'graph';
end

if nargin < 3
    options = [];
end

if nargin<4
    extstruct = yalmip('extstruct',extvar);
elseif isempty(extstruct)
    extstruct = yalmip('extstruct',extvar);
end

if isempty(extstruct)
    error('This is not a nonlinear operator variable');
end

fcn = extstruct.fcn;
try
    [F,properties,arguments] = feval(fcn,method,extstruct.var,extstruct.arg{1:end-1});
catch
    error(['Failed when trying to create a model for the "' extstruct.fcn '" operator']);
end

% These fields are not official, or not required yet
if ~isempty(properties)
    properties.name = fcn;
    if ~isfield(properties,'derivative')
        properties.derivative = [];
    end
    if ~isfield(properties,'models')
        properties.models = getvariables(extstruct.var);
    end
    if ~isfield(properties,'convexhull')
        properties.convexhull = [];
    end
    if ~isfield(properties,'bounds')
        properties.bounds = [];
    end
    if ~isfield(properties,'domain')
        properties.domain = [-inf inf];
    end
    if ~isfield(properties,'range')
        if strcmpi(properties.definiteness,'positive')
            properties.range = [0 inf];
        elseif strcmpi(properties.definiteness,'negative')
            properties.range = [-inf 0];
        else
            properties.range = [-inf inf];
        end
    end
    if ~isfield(properties,'model')
        properties.model = 'unspecified';
    end
end

% Normalize the callback expression and check for some obsoleted stuff
if ~isempty(properties)
    if isequal(properties.model,'callback')
        F_normalizing = NormalizeCallback(method,extstruct.var,extstruct.arg{:});
        F = F + F_normalizing;
    end
    if length(extstruct.computes)>1
        properties.models = extstruct.computes;
    end
    if ~any(strcmpi(properties.convexity,{'convex','concave','none'}))
        disp('More cleaning, strange convextiy returned...Report bug')
        error('More cleaning, strange convextiy returned...Report bug')
    end
end

% This is useful in MPT
if ~isempty(F)
    F = tag(F,['Expansion of ' extstruct.fcn]);
end
