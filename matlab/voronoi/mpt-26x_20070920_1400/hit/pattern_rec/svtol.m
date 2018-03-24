function tol = svtol(C)
%SVTOL Tolerance for Support Vectors
%
%  Usage: tol = svtol(C)
%
%  Parameters: C      - upper bound 
%
%  Author: Steve Gunn (srg@ecs.soton.ac.uk)


  if (nargin ~= 1) % check correct number of arguments
    help svtol
  else

    % tolerance for Support Vector Detection
    if C==Inf 
      tol = 1e-5;
    else
      tol = C*1e-6;
    end
    
  end
