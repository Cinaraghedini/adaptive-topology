function varargout = sumabsk(varargin)
% SUMABSK  Returns sum of k largest (by magnitude) (eigen-)values.
%
% s = SUMABSK(X,k) 
%
% For a vector X, SUMABSK returns the sum of the k largest (by magnitude) elements.
%
% For a symmetric matrix X, SUMABSK returns the sum of the k largest (by magnitude) eigen-values.
%
% See also SUMK

% Author Johan Löfberg
% $Id: sumabsk.m,v 1.3 2007/08/02 19:17:36 joloef Exp $

% ***************************************************
% This file defines a nonlinear operator for YALMIP
%
% It can take three different inputs
% For double inputs, it returns standard double values
% For sdpvar inputs, it genreates a an internal variable
% When first input is 'model' it generates the epigraph
%
% ***************************************************
switch class(varargin{1})

    case 'double' % What is the numerical value of this argument (needed for displays etc)
        if nargin == 1
            error('sumabsk needs two arguments');
        else
            X = varargin{1};
            [n,m] = size(X);
            if (min(n,m) > 1 & ~issymmetric(X))
                error('sumabsk can only be applied on vectors and symmetric matrices');
            else            
                k = min(length(X),varargin{2});
                if min(n,m)==1                    
                    sorted = sort(X);                   
                else
                    sorted = sort(eig(X));
                end
                varargout{1} = sum(sorted(max(1,end-k+1):end));
            end
        end

    case 'sdpvar' % Overloaded operator for SDPVAR objects. Pass on args and save them.
        X = varargin{1};
        [n,m] = size(X);
        if (min(n,m) > 1 & ~issymmetric(X))
             error('sumabsk can only be applied on vectors and symmetric matrices');
        else
            if nargin < 2
                error('sumabsk needs two arguments');
            else
                varargout{1} = yalmip('define',mfilename,varargin{:});
            end
        end

    case 'char' % YALMIP sends 'model' when it wants the epigraph or hypograph
        if isequal(varargin{1},'graph')
            t = varargin{2}; % Second arg is the extended operator variable
            X = varargin{3}; % Third arg and above are the args user used when defining t.
            k = min(varargin{4},length(X));
            [n,m] = size(X);
            Z = sdpvar(n,m);
            s = sdpvar(1,1);
            if min(n,m)==1
                varargout{1} = set(t-k*s-sum(Z) > 0) + set(Z > 0) + set(Z+s > X > -Z-s);
            else
                varargout{1} = set(t-k*s-trace(Z) > 0) + set(Z > 0) + set(Z+s*eye(n) > X > -Z-s*eye(n));
            end
            varargout{2} = struct('convexity','convex','monotonicity','none','definiteness','none','model','graph');
            varargout{3} = X;
        else
        end
    otherwise
end
