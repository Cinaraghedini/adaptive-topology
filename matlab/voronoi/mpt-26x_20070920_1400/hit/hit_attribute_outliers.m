%HIT_ATTRIBUTE_OUTLIERS Attribute outliers to the modes and update inliers accordingly
% -------------------------------------------------------------------------
% DESCRIPTION
% -------------------------------------------------------------------------
%
% This is a batch file.
% 
% -------------------------------------------------------------------------
% INPUT
% -------------------------------------------------------------------------
% (none)
%
% -------------------------------------------------------------------------
% OUTPUT                                                                                                   
% -------------------------------------------------------------------------
% (none)

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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INSERT HERE YOUR PIECE OF CODE FOR ATTRIBUTING
% OUTLIERS TO MODE DATASETS (AND UPDATE INLIERS ACCORDINGLY)
% At present, outliers are attributed to the most likely mode
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%% Now attribute outliers, if any to the most_likely mode ...
% i.e. the one that minimizes the squared distance
% a possible pitfall happens when there is the same PV for different modes.
% In this case, the outlier may be attributed to the wrong mode
% since information about the "relative position" of modes and 
% outliers is not taken into account in the attribution rule.

% % outl_attributed: flag describing if some outliers have been attributed or not
outl_attributed=0;
out_pos=setdiff([1:ndata],inliers.pos);
if ~isempty(out_pos)
    outl_attributed=1;
    noutl=length(out_pos);
    for kk=1:noutl
        ytempp=yid(out_pos(kk));
        Xtempp=Xid(out_pos(kk),:);
        distances=[];
        for ii=1:s
            distances=[distances;(ytempp-[Xtempp , 1]*idmodes1.par{ii}).^2];
        end %ends for ii
        [tempp,newi]=min(distances); %newi store the index of the most likely mode
        % add the point to the exchange list
        inliers.pos=[inliers.pos;out_pos(kk)];
        inliers.class=[inliers.class;newi];
    end
else
    fprintf('\n No datapoint has been reattributed\n');
end
