function M = incidence_matrix(g)
% incidence_matrix(g) --- return the vertex/edge incidence matrix.
% We return M an nv-by-ne matrix whose ij entry is 1 if vertex i is an
% end point of edge j.

[n,m] = size(g);
e = edges(g);

M = logical(zeros(n,m));

for k=1:m
    a = e(k,1);
    b = e(k,2);
    M(a,k)=1;
    M(b,k)=1;
end