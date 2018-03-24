function [y,region]=hit_pwa(theta,regions,point)
%HIT_PWA Evaluate a PWA function.
%
% -------------------------------------------------------------------------
% DESCRIPTION
% -------------------------------------------------------------------------
%
% [y,region]=hit_pwa(theta,regions,point)
% 
% -------------------------------------------------------------------------
% INPUT
% -------------------------------------------------------------------------
% theta{i}: PV of the i-th mode.
%
% regions(i): polyhedron describing the region of the i-th mode. 
%
% point: vector describing the point at which the map will be evaluated (it
% can be either a row or a column).
% 
% -------------------------------------------------------------------------
% OUTPUT                                                                                                   
% -------------------------------------------------------------------------
% y: value of the PWA map. It is NaN if the point does not fall in at least
% one region. This allows one to skip the point when using plot or plot3.
%
% region: number of the region the point belongs to. It is NaN if the point
% does not fall in at least one region.

% Copyright is with the following author:
%
% (C) 2005 Giancarlo Ferrari Trecate,
%         giancarlo.ferrari@unipv.it
% -------------------------------------------------------------------------
% Legal note:
%     This program is free software; you can redistribute it and/or
%     modify it under the terms of the GNU General Public
%     License as published by the Free Software Foundation; either
%     version 2.1 of the License, or (at your option) any later version.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%     General Public License for more details.
%
%     You should have received a copy of the GNU General Public
%     License along with this library; if not, write to the
%     Free Software Foundation, Inc.,
%     59 Temple Place, Suite 330,
%     Boston, MA  02111-1307  USA
%
%
% -------------------------------------------------------------------------

nreg=length(theta); % # of regions and models

% Return NaN if the x-point do not fall in at least one region
% this correspond to avoid plotting the point when using
% the commands plot or plot3 !

% Note: if a point is on the boundary of two regions,
% say region(k) and region (j) it is assumed that the
% the mode with lower index is active.
% For instance, if j<k, mode(j) is active

point=point(:);
y=NaN;
region=NaN;
k=1;
while isnan(region) & k<=nreg
    test=isinside(regions(k),point); % 1 if point belongs to the kth region
    if test==1
        region=k;
        th=theta{k};
        th=th(:);
        y=th'*[point;1];  % compute the pwa function
    end
    k=k+1;
end
