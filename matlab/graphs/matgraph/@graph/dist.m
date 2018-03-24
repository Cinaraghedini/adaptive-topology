function d = dist(g,v,w)
% dist(g,v,w) and dist(g,v) --- find distance(s) between vertices
% The form dist(g,v,w) finds the distance between vertices v and w.
% The form dist(g,v) returns a vector distance from v to all other vertices
% in the graph.
% The form dist(g) returns an n-by-n matrix whose ij entry is the distance
% between vertices i and j.
%
% The code for all pairs shortest distance was written by Michael Kleder
% and found on the Mathwork's MATLAB Central file exchange.


if nargin == 3
    p = find_path(g,v,w);
    if isempty(p)
        d = inf;
    else
        d = length(p)-1;
    end
    return
end

n = nv(g);


if nargin==1
    d = allspath(matrix(g));
    return
end



d = ones(1,n) * Inf;

if (v<0) | (v>n)
    return
end

d(v) = 0;
d = d';

t = graph;
bfstree(t,g,v);
sparse(t);

vec = zeros(n,1);
vec(v) = 1;
A = double(matrix(t));
vec_sum = sum(vec);

for k=1:n
    vec2 = A*vec > 0;
    vec2 = vec2 & (d==Inf);
    d(vec2) = k;
    vec = vec2 + vec > 0;
    new_vec_sum = sum(vec);
    if (new_vec_sum == vec_sum)
        break;
    else
        vec_sum = new_vec_sum;
    end
end
d = d';
free(t);    


% this subroutine was written by Michael Kleder and was found on The
% Mathwork's MATLAB Central file exchange.

function B = allspath(A)
B=full(double(A));
B(B==0)=Inf;
C=ones(size(B));
iter=0;
while any(C(:))
    C=B;
    B=min(B,squeeze(min(repmat(B,[1 1 length(B)])+...
        repmat(permute(B,[1 3 2]),[1 length(B) 1]),[],1)));
    C=B-C;
end
B(logical(eye(length(B))))=0;
return