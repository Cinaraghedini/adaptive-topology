function x = uncertain(x)
%UNCERTAIN Declares a variable as uncertain
%
%   F = UNCERTAIN(x) is used to describe the set of uncertain variables
%   in an uncertain program, as an alternative to using the 4th input in
%   SOLVEROBUST.
%
%   If an uncertain multi-parametric problem is solved, UNCERTAIN has to be
%   used to declare the set of uncertain variables (since the fourt
%   argument in solvemp is devoted to the paramtric variable)
%
%   INPUT
%    x : SDPVAR object
%
%   OUTPUT
%    F : SET object
%
%   See also SOLVEROBUST, ROBUSTIFY

% Author Johan Löfberg
% $Id: uncertain.m,v 1.3 2006/08/18 15:01:04 joloef Exp $

x.typeflag = 15;
x = lmi(x);