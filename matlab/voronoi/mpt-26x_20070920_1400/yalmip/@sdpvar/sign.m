function varargout = SIGN(varargin)
%SIGN (overloaded)

% Author Johan Löfberg
% $Id: sign.m,v 1.3 2007/07/30 15:12:14 joloef Exp $
switch class(varargin{1})

    case 'double' % What is the numerical value of this argument (needed for displays etc)
        % SHOULD NEVER HAPPEN, THIS SHOULD BE CAUGHT BY BUILT-IN
        error('Overloaded SDPVAR/NORM CALLED WITH DOUBLE. Report error')

    case 'sdpvar' % Overloaded operator for SDPVAR objects. Pass on args and save them.
        if length(varargin{1}) == 1
            varargout{1} = yalmip('addextendedvariable',mfilename,varargin{1});
        else
            y = [];
            n = size(varargin{1},1);
            m = size(varargin{1},2);
            varargin{1} = reshape(varargin{1},n*m,1);
            for i = 1:prod(size(varargin{1}))
                inparg = extsubsref(varargin{1},i);
                if isa(inparg,'double')
                    y = sign(inparg);
                else
                    y = [y yalmip('addextendedvariable',mfilename,inparg)];
                end
            end
            y = reshape(y,n,m);
            varargout{1} = y;
        end

    case 'char' % YALMIP sends 'model' when it wants the epigraph or hypograph
        switch varargin{1}
            case {'graph','milp'}

                t = varargin{2};
                X = varargin{3};

                d = binvar(1,1);
                [M,m] = derivebounds(X);
                F = [X >= d*m,X <=(1-d)*M, t == 1-2*d];

                varargout{1} = F;
                varargout{2} = struct('convexity','none','monotonicity','none','definiteness','none');
                varargout{3} = X;

            otherwise
                error('SDPVAR/SIGN called with CHAR argument?');
        end
    otherwise
        error('Strange type on first argument in SDPVAR/SIGN');
end
