function varargout = acos(varargin)
%ACOS (overloaded)

% Author Johan Löfberg
% $Id: acos.m,v 1.4 2007/08/02 18:16:26 joloef Exp $
switch class(varargin{1})

    case 'double'
        error('Overloaded SDPVAR/ACOS CALLED WITH DOUBLE. Report error')

    case 'sdpvar'
        varargout{1} = InstantiateElementWise(mfilename,varargin{:});

    case 'char'

        operator = struct('convexity','none','monotonicity','decreasing','definiteness','none','model','callback');
        operator.convexhull = [];
        operator.bounds = @bounds;
        operator.derivative = @(x)(-(1 - x.^2).^-0.5;
        operator.range = [-pi/2 pi/2];
        operator.domain = [-1 1];

        varargout{1} = [];
        varargout{2} = operator;
        varargout{3} = varargin{3};

    otherwise
        error('SDPVAR/ACOS called with CHAR argument?');
end

function [L,U] = bounds(xL,xU)
L = acos(xU);
U = acos(xL);