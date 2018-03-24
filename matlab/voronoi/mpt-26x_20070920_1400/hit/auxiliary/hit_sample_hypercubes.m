function X=hit_sample_hypercubes(S,n,Bigregion)
%HIT_SAMPLE_HYPERCUBES Generate random points in a list of n-dimensional hypercubes intersected with a bounding polytope.
%
% -------------------------------------------------------------------------
% DESCRIPTION
% -------------------------------------------------------------------------
%
% X=hit_sample_hypercubes(S,n,Bigregion)
%
% Example: X=sampling_hypercubes...
% ({{1,[0 0]',20},{2,[1,-2]',500}},2,Bigregion)
% 
% -------------------------------------------------------------------------
% INPUT
% -------------------------------------------------------------------------
% S{i} is a cell {l,c,np}
%   l: lenght on each edge of the i-th hypercube.
%   c: center of the i-th hypercube.
%   np: number of points to be extracted in the i-th hypercube.
%
% n: dimension of each hypercube.
%
% Bigregion: n-dimensional bounding polytope.
% 
% -------------------------------------------------------------------------
% OUTPUT                                                                                                   
% -------------------------------------------------------------------------
% X(i,:) i-th sampled point.
%
% -------------------------------------------------------------------------
% COMMENTS                                                                                              
% -------------------------------------------------------------------------
% This function does not work if the bounding polytope and *some* hypercube
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

Ncubes=length(S);
X=[];
eps1=1e-6;
A=[eye(n);-eye(n)];
for i=1:Ncubes
    % check that Bigregion and the i-th hypercube
    % HAVE nonempty intersection
    tl=S{i}{1};
    tc=S{i}{2};
    tc=tc(:);
    tnp=S{i}{3};
    minimals=tc-tl/2; minimals=minimals(:);
    maximals=tc+tl/2; maximals=maximals(:);
    I=polytope(A,[maximals;-minimals]);
    P=Bigregion&I;
    [xCheb,RCheb] = chebyball(P);
    if RCheb<1e-6
        fprintf('\n ======================================================== ');
        fprintf('\n The intersection of the %d-th interval with the Bigregion',i);
        fprintf('\n is too small (Chebyshev radius <%f', eps1);
        fprintf('\n ======================================================== \n');
        error('hit_sampling_boxes:end','hit_sampling_boxes ends here.')
    end
    j=0;
    while(j<tnp)
        point=(rand(n,1)-0.5)*tl+tc;
        isin= isinside(P,point);
        if isin
            j=j+1;
            X=[X;point'];
        end
    end
end

    