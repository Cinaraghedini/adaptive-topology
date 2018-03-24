function issym=issymmetric(X)
%ISSYMMETRIC Check if variable is symmetric

% Author Johan Löfberg 
% $Id: issymmetric.m,v 1.2 2004/07/02 08:17:31 johanl Exp $

[n,m] = size(X);
issym = 0;

if (n==m)
  issym = norm(X-X',1)<1e-10; % Should be scaled with size maybe
end

  
  
