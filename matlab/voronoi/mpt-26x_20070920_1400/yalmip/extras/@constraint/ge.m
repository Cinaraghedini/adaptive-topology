function F = gt(X,Y)
% Internal class for constraint lists

% Author Johan Löfberg
% $Id: ge.m,v 1.1 2004/06/17 08:40:02 johanl Exp $

superiorto('sdpvar');
superiorto('double');

try
    % Try to evaluate
    if isa(X,'constraint')
        % (z > w) > y
        Z = X.List{end} - Y;
        F = X;
        F.List{end+1} = '>';
        F.List{end+1} = Y;
        F.Evaluated{end+1} = Z;
        F.strict(end+1) = 0;
    else
        % z > (w > y)
        Z = X - Y.List{1};
        F = Y;
        F.List = {X,'>',F.List{:}};
        F.Evaluated = {Z,F.Evaluated{:}};
        F.strict = [1 F.strict];               
    end
catch
    error(lasterr);
end

