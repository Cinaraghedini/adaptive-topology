function k = svkernel(ker,u,v)
%SVKERNEL Kernel for Support Vector Methods
%
%  Usage: k = svkernel(ker,u,v)
%
%  Parameters: ker - kernel type
%              u,v - kernel arguments
%
%  Values for ker: 'linear'  -
%                  'poly'    - p1 is degree of polynomial
%                  'rbf'     - p1 is width of rbfs (sigma)
%                  'sigmoid' - p1 is scale, p2 is offset
%                  'spline'  -
%                  'bspline' - p1 is degree of bspline
%                  'fourier' - p1 is degree
%                  'erfb'    - p1 is width of rbfs (sigma)
%                  'anova'   - p1 is max order of terms
%              
%  Author: Steve Gunn (srg@ecs.soton.ac.uk)

  if (nargin < 1) % check correct number of arguments
     help svkernel
  else
     
    global p1 p2;

    % could check for correct number of args in here
    % but will slow things down further
    switch lower(ker)
      case 'linear'
        k = u*v';
      case 'poly'
        k = (u*v' + 1)^p1;
      case 'rbf'
        k = exp(-(u-v)*(u-v)'/(2*p1^2));
      case 'erbf'
        k = exp(-sqrt((u-v)*(u-v)')/(2*p1^2));
      case 'sigmoid'
        k = tanh(p1*u*v'/length(u) + p2);
      case 'fourier'
        z = sin(p1 + 1/2)*2*ones(length(u),1);
        i = find(u-v);
        z(i) = sin(p1 + 1/2)*(u(i)-v(i))./sin((u(i)-v(i))/2);
        k = prod(z);
      case 'spline'
        z = 1 + u.*v + u.*v.*min(u,v) - ((u+v)/2).*(min(u,v)).^2 + (1/3)*(min(u,v)).^3;
        k = prod(z);
      case {'curvspline','anova'}
        z = 1 + u.*v + (1/2)*u.*v.*min(u,v) - (1/6)*(min(u,v)).^3;
        k = prod(z);

% - sum(u.*v) - 1; 
%        z = 1 + u.*v + (1/2)*u.*v.*min(u,v) - (1/6)*(min(u,v)).^3;
%        k = prod(z); 
%        z = (1/2)*u.*v.*min(u,v) - (1/6)*(min(u,v)).^3;
%        k = prod(z); 

      case 'bspline'
        z = 0;
        for r = 0: 2*(p1+1)
          z = z + (-1)^r*binomial(2*(p1+1),r)*(max(0,u-v + p1+1 - r)).^(2*p1 + 1);
        end
        k = prod(z);
      case 'anovaspline1'
        z = 1 + u.*v + u.*v.*min(u,v) - ((u+v)/2).*(min(u,v)).^2 + (1/3)*(min(u,v)).^3;
        k = prod(z); 
      case 'anovaspline2'
        z = 1 + u.*v + (u.*v).^2 + (u.*v).^2.*min(u,v) - u.*v.*(u+v).*(min(u,v)).^2 + (1/3)*(u.^2 + 4*u.*v + v.^2).*(min(u,v)).^3 - (1/2)*(u+v).*(min(u,v)).^4 + (1/5)*(min(u,v)).^5;
        k = prod(z);
      case 'anovaspline3'
        z = 1 + u.*v + (u.*v).^2 + (u.*v).^3 + (u.*v).^3.*min(u,v) - (3/2)*(u.*v).^2.*(u+v).*(min(u,v)).^2 + u.*v.*(u.^2 + 3*u.*v + v.^2).*(min(u,v)).^3 - (1/4)*(u.^3 + 9*u.^2.*v + 9*u.*v.^2 + v.^3).*(min(u,v)).^4 + (3/5)*(u.^2 + 3*u.*v + v.^2).*(min(u,v)).^5 - (1/2)*(u+v).*(min(u,v)).^6 + (1/7)*(min(u,v)).^7;
        k = prod(z);
      case 'anovabspline'
        z = 0;
        for r = 0: 2*(p1+1)
          z = z + (-1)^r*binomial(2*(p1+1),r)*(max(0,u-v + p1+1 - r)).^(2*p1 + 1);
        end
        k = prod(1 + z);
      otherwise
        k = u*v';
    end

  end
