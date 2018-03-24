function [xunit,yunit] = GeneratePointsCircle(radius, cx, cy)
%function [xunit,yunit] = GeneratePointsCircle(radius, cx, cy)
% Generate points that in circle (for further plot)
% Parameters: radius: degrees (float)
%             cx: degrees (float)
%             cy: degrees (float)
% Output: xunit = xRange
%         yunut = yRange
% Technological Institute of Aeronautics
% Author: Nicolas Pereira Borges - nicolas@ita.br
% Date: 20/09/2016

   th = 0:pi/50:2*pi;
   xunit = radius * cos(th) + cx;
   yunit = radius * sin(th) + cy;
   

end