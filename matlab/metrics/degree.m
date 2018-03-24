%  degree.m
%  Ray  Byrne
%  August  2,  2005
%  degree.m returns the  D  matrix where  the  (i,i) term  is
%  the  degree of  the  ith  node.   The  input  is  the  adjacency matrix A
function D  =  degree(A)
degrees =  sum(A);
N  =  length(degrees);
D  =  eye(N);
for  i=1:N
    D(i,i) =  degrees(i);
end  %  for
return