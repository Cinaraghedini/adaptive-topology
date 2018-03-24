% X = cell2diag(x)
% Creates diag(x{1},x{2},..x{length(x)}) block diagonal matrix.
%
% SEE ALSO cellK
function X=cell2diag(x)
% (C) Jos F. Sturm, 2001.
 sizex = zeros(1+length(x),2);
 for k=1:length(x)
     sizex(k+1,:) = size(x{k});
 end
 is = cumsum(sizex(:,1));
 js = cumsum(sizex(:,2));
 X = sparse([],[],[],is(end),js(end));
 for k=1:length(x)
     X(is(k)+1:is(k+1),js(k)+1:js(k+1)) = x{k};
 end
