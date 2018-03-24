function varargout = log(varargin)
%LOG (overloaded)

% Author Johan Löfberg
% $Id: log.m,v 1.16 2007/08/17 19:11:42 joloef Exp $
switch class(varargin{1})

    case 'double'
        error('Overloaded SDPVAR/LOG CALLED WITH DOUBLE. Report error')

    case 'sdpvar'
        varargout{1} = InstantiateElementWise(mfilename,varargin{:});

    case 'char'

        X = varargin{3};      
        F = set(X > 1e-8);

        operator = struct('convexity','concave','monotonicity','increasing','definiteness','none','model','callback');
        operator.convexhull = @convexhull;
        operator.bounds = @bounds;
        operator.derivative = @(x)(1./(abs(x)+eps));

        varargout{1} = F;
        varargout{2} = operator;
        varargout{3} = X;

    otherwise
        error('SDPVAR/LOG called with CHAR argument?');
end

function [L,U] = bounds(xL,xU)
if xL <= 0
    % The variable is not bounded enough yet
    L = -inf;
else
    L = log(xL);
end
if xU < 0
    % This is an infeasible problem
    L = inf;
    U = -inf;
else
    U = log(xU);
end

function [Ax, Ay, b] = convexhull(xL,xU)
if xL <= 0
    fL = inf;
else
    fL = log(xL);
end
fU = log(xU);
dfL = 1/(xL);
dfU = 1/(xU);
xM = (xU - xL)/(fU-fL);
fM = log(xM);
dfM = 1/xM;

[Ax,Ay,b] = convexhullConcave(xL,xM,xU,fL,fM,fU,dfL,dfM,dfU);
remove = isinf(b) | isinf(Ax) | isnan(b);
if any(remove)
    remove = find(remove);
    Ax(remove)=[];
    b(remove)=[];
    Ay(remove)=[];
end