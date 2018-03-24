function X=hit_sample_intervals(S, Bigregion)
%HIT_SAMPLE_INTERVALS Generate random points in a list of intervals intersected with a bounding interval.
%
% -------------------------------------------------------------------------
% DESCRIPTION
% -------------------------------------------------------------------------
%
% X=hit_sample_intervals(S, Bigregion)
%
% Example: X=sample_intervals...
% ({{[1,2],20},{[-2,0],500}},Bigregion)
% 
% -------------------------------------------------------------------------
% INPUT
% -------------------------------------------------------------------------
% S{i} is a cell {int,np}
%   int: [a,b] i-th interval.
%   np: number of points to be extracted in the i-th interval (intersected
%   with Bigregion).
%
% Bigregion: bounding interval.
% 
% -------------------------------------------------------------------------
% OUTPUT                                                                                                   
% -------------------------------------------------------------------------
% X(i) i-th sampled point.
%
% -------------------------------------------------------------------------
% COMMENTS                                                                                              
% -------------------------------------------------------------------------
% This function does not work if the bounding interval and *some* interval
% have an empty intersection. 

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

Nbox=length(S);
X=[];
eps1=1e-6;
for i=1:Nbox
    % check that Bigregion and the i-th interval
    % HAVE nonempty intersection
    I=polytope(S{i}{1}(:));
    P=Bigregion&I;
    [xCheb,RCheb] = chebyball(P);
    if RCheb<1e-6
        fprintf('\n ======================================================== ');
        fprintf('\n The intersection of the %d-th interval with the Bigregion',i);
        fprintf('\n is too small (Chebyshev radius <%f', eps1);
        fprintf('\n ======================================================== \n');
        error('hit_sampling_intervals:end','hit_sampling_intervals ends here.')
    end
    tnp=S{i}{2};
    j=0;
    %center of the interval
    V=extreme(P);
    V=sort(V);
    tc=[(V(2)-V(1))/2+V(1)];
    tl=[(V(2)-V(1))]; %length of each interval
    for j=1:tnp
        point=(rand(1,1)-0.5).*tl+tc;
        X=[X;point];
    end
end

    