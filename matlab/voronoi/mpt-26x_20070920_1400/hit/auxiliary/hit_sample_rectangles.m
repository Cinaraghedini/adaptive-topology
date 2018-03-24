function X=hit_sample_rectangles(S, Bigregion)
%HIT_SAMPLE_RECTANGLES Generate random points in a list of 2D rectangles intersected with a bounding polytope.
%
% -------------------------------------------------------------------------
% DESCRIPTION
% -------------------------------------------------------------------------
%
% X=hit_sample_rectangles(S, Bigregion)
%
% Example: X=hit_sample_rectangles...
% ({{[1,2],[0,1],20},{[2,1],[-2,0],500}},Bigregion)
% 
% -------------------------------------------------------------------------
% INPUT
% -------------------------------------------------------------------------
% S{i} is a cell {int1,int2,np}
%   int1: [a,b] interval on the first coordinate describing the i-th
%   rectangle.
%   int2: [a',b'] interval on the second coordinate describing the i-th
%   rectangle.
%   np: number of points to be extracted in the i-th rectangle (intersected
%   with Bigregion).
%
% Bigregion: bounding polytope.
% 
% -------------------------------------------------------------------------
% OUTPUT                                                                                                   
% -------------------------------------------------------------------------
% X(i,:) i-th sampled point.
%
% -------------------------------------------------------------------------
% COMMENTS                                                                                              
% -------------------------------------------------------------------------
% This function does not work if the bounding polytope and *some* rectangle
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

A=[eye(2);-eye(2)];
for i=1:Nbox
    % check that Bigregion and the i-th interval
    % HAVE nonempty intersection
    I1=sort(S{i}{1});
    I2=sort(S{i}{2});
    minimals=[I1(1),I2(1)]';
    maximals=[I1(2),I2(2)]';
    I=polytope(A,[maximals;-minimals]);
    P=Bigregion&I;
    [xCheb,RCheb] = chebyball(P);
    if RCheb<1e-6
        fprintf('\n ======================================================== ');
        fprintf('\n The intersection of the %d-th interval with Bigregion',i);
        fprintf('\n is too small (Chebyshev radius <%f', eps1);
        fprintf('\n ======================================================== \n');
        error('hit_sampling_boxes:end','hit_sampling_boxes ends here.')
    end
    tnp=S{i}{3};
    j=0;
    tc=[(I1(2)-I1(1))/2+I1(1),(I2(2)-I2(1))/2+I2(1)]'; %center of the box
    tl=[(I1(2)-I1(1)),(I2(2)-I2(1))]'; %length of each edge
    while(j<tnp)
        % generate a point in the box
        point=(rand(2,1)-0.5).*tl+tc;
        % if the point belongs also to Bigregion, keep it 
         isin= isinside(P,point);
        if isin
            j=j+1;
        X=[X;point'];
        end
    end
end

