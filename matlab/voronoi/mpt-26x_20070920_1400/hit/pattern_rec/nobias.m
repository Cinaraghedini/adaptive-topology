function nb = nobias(ker)
%NOBIAS Returns true if SVM kernel has no implicit bias
%
%  Usage: nb = nobias(ker)
%
%  Parameters: ker - kernel type
%              
%  Author: Steve Gunn (srg@ecs.soton.ac.uk)

  if (nargin ~= 1) % check correct number of arguments
     help nobias
  else
     
    switch lower(ker)
      case {'linear','sigmoid'}
%,'anovaspline1','anovaspline2','anovaspline3'}
        nb = 1;
      otherwise
        nb = 0;
    end

  end
