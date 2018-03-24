function X = not(X)
% Internal class for constraint list

% Author Johan Löfberg
% $Id: not.m,v 1.1 2004/06/17 08:40:02 johanl Exp $

superiorto('sdpvar');
superiorto('double');

% Try to evaluate

switch X.List{2}
    case '<'
        X.List{2} = '>';
    case '>'
        X.List{2} = '<';
    otherwise
        error
end

X.Evaluated{1} = -X.Evaluated{1};