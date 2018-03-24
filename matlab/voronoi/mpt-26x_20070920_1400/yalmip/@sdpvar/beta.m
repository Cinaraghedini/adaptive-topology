function varargout=beta(varargin)
%BETA (overloaded)

% Author Johan Löfberg
% $Id: beta.m,v 1.3 2007/08/02 18:16:26 joloef Exp $
switch class(varargin{1})

    case 'sdpvar'
        varargout{1} = InstantiateElementWise(mfilename,varargin{:});

    case 'char'

        X = varargin{3};
        F = set(X > eps);
        operator = struct('convexity','convex','monotonicity','decreasing','definiteness','positive');
        varargout{1} = F;
        varargout{2} = operator;
        varargout{3} = varargin{3};

    otherwise
        error('Strange type on first argument in SDPVAR/BETA');
end