function varargout = plot(varargin)
% Internal class for constraint list

% Author Johan Löfberg
% $Id: plot.m,v 1.1 2007/02/28 16:20:33 joloef Exp $

varargin{1} = lmi(varargin{1});
plot(varargin{:});
