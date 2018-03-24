function sys = dual(X)
%DUAL Extract dual variable
%   
%   Z = DUAL(F)     Returns the dual variable for the constraint F
% 
%   See also   SET, SOLVESDP
  
% Author Johan Löfberg
% $Id: dual.m,v 1.3 2005/02/04 10:10:26 johanl Exp $

nlmi = size(X.clauses,2);

% Is it an empty SET
if (nlmi == 0) 
    sys = [];
    return
end

if nlmi>1
    error('Dual not applicable on list of constraints. Use dual(F(index)) or dual(F(''tag''))')
end

% Get dual from repospitory
sys = yalmip('dual',X.LMIid);

% If no dual available, returns NaNs with correct dimension
if isempty(sys)
    sys = real(double(X))+NaN;
end