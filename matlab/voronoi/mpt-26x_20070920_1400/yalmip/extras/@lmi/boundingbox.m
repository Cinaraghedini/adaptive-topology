function [L,U,x] = boundingbox(F,ops)
%BOUINDINGBOX Computes bounding box of a SET object
%
% [L,U,x] = boundingbox(F)

% Author Johan Löfberg
% $Id: boundingbox.m,v 1.1 2004/12/08 00:07:15 johanl Exp $

x = recover(depends(F));

if nargin < 2
    ops = sdpsettings('verbose',0);    
end

for i = 1:length(x);
    sol = solvesdp(F,x(i),ops);
    L(i,1) = double(x(i));
    sol = solvesdp(F,-x(i),ops);
    U(i,1) = double(x(i));
end