function varargout = cos(varargin)
%COS (overloaded)

% Author Johan Löfberg
% $Id: cos.m,v 1.12 2007/08/29 19:05:54 joloef Exp $
switch class(varargin{1})

    case 'double'
        error('Overloaded SDPVAR/COS CALLED WITH DOUBLE. Report error')

    case 'sdpvar'
        varargout{1} = InstantiateElementWise(mfilename,varargin{:});
        
    case 'char'
        
        operator = struct('convexity','none','monotonicity','none','definiteness','none','model','callback');
        operator.bounds     = @bounds;
        operator.derivative = @(x)(-sin(x));
        operator.range = [-1 1];
        operator.convexhull = @convexhull;
 
        varargout{1} = [];
        varargout{2} = operator;
        varargout{3} = varargin{3};

    otherwise
        error('SDPVAR/COS called with CHAR argument?');
end

function [L,U] = bounds(xL,xU)
if xU-xL >= 2*pi
    L = -1;
    U = 1;
else
    xL = xL + pi/2;
    xU = xU + pi/2;
    n = floor(( (xL + xU)/2/(2*pi)));
    xL = xL - n*2*pi;
    xU = xU - n*2*pi;
    yL = sin(xL);
    yU = sin(xU);
    L = min([yL yU]);
    U = max([yL yU]);
    if (xL<pi/2 & xU>pi/2) |  (xL<-3*pi/2 & xU>-3*pi/2)
        U = 1;
    end
    if (xL < 3*pi/2 & xU > 3*pi/2) | (xL < -pi/2 & xU > -pi/2)
        L = -1;
    end
end

function [Ax, Ay, b] = convexhull(xL,xU)
if cos(xL)>=0 & cos(xU)>=0 & xU-xL<pi
    fL = cos(xL);
    fU = cos(xU);
    dfL = -sin(xL);
    dfU = -sin(xU);
    [Ax,Ay,b] = convexhullConcave(xL,xU,fL,fU,dfL,dfU);
elseif cos(xL)<=0 & cos(xU)<=0 & xU-xL<pi
    fL = cos(xL);
    fU = cos(xU);
    dfL = -sin(xL);
    dfU = -sin(xU);
    [Ax,Ay,b] = convexhullConvex(xL,xU,fL,fU,dfL,dfU);
else
    z = linspace(xL,xU,100);
    fz = cos(z);
    [minval,minpos] = min(fz);
    [maxval,maxpos] = max(fz);
    xtestmin = linspace(z(max([1 minpos-5])),z(min([100 minpos+5])),100);
    xtestmax = linspace(z(max([1 maxpos-5])),z(min([100 maxpos+5])),100);

    fz1 = cos(xtestmin);
    fz2 = cos(xtestmax);
    z = [z(:);xtestmin(:);xtestmax(:)];
    fz = [fz(:);fz1(:);fz2(:)];
    [z,sorter] = sort(z);
    fz = fz(sorter);
    [z,ii,jj]=unique(z);
    fz = fz(ii);

    k1 = max((fz(2:end)-fz(1))./(z(2:end)-xL))+1e-12;
    k2 = min((fz(2:end)-fz(1))./(z(2:end)-xL))-1e-12;
    k3 = min((fz(1:end-1)-fz(end))./(z(1:end-1)-xU))+1e-12;
    k4 = max((fz(1:end-1)-fz(end))./(z(1:end-1)-xU))-1e-12;
    Ax = [-k1;k2;-k3;k4];
    Ay = [1;-1;1;-1];
    b =  [k1*(-z(1)) + fz(1);-(k2*(-z(1)) + fz(1));k3*(-z(end)) + fz(end);-(k4*(-z(end)) + fz(end))];
end