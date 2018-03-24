function [inliers_new,s_new]=hit_check_inliers(inliers,s,ndim)
%HIT_CHECK_INLIERS Remove modes, if necessary, from a list of inliers
%
% -------------------------------------------------------------------------
% DESCRIPTION
% -------------------------------------------------------------------------
%
% [inliers_new,s_new]=hit_check_inliers(inliers,s,ndim)
% 
% -------------------------------------------------------------------------
% INPUT
% -------------------------------------------------------------------------
% inliers: structure containing information on the inliers
%
%   inliers.pos(j): index of the j-th inlier in Xid and yid
%	(i.e.Xid(inliers.pos(j),:) and yid(inliers.pos(j)) are the j-th
%	inliers).
%   inliers.class(j): classification of the j-th inlier.
%
% s: maximal number of modes the inliers are attributed to.
%
% ndim: dimension of each inlier.
%
% -------------------------------------------------------------------------
% OUTPUT                                                                                                   
% -------------------------------------------------------------------------
% inliers_new: structure containing information on the inliers. It has the same
% fields of inliers.
%
% s_new: updated number of modes.

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


global idpar ;
if ~isstruct(idpar)
    hit_error;
end


% Remove clusters that have less than idpar.c points and decrement
% idpar.s
% The points in such clusters are removed from the list of inliers 

% mod_cut stores how many clusters are removed
mod_cut=0;
ind_inl_new=[];
ind_class_inl_new=[];
for i=1:s
    pos=find(inliers.class==i);
    %check if there are less points than n. of parameters
    if length(pos)<idpar.c 
        % cut the points
        fprintf('\n WARNING: Mode %d has less than %d points (the size of an LD) \n',i,idpar.c )
        fprintf('Its points are discarded and the number of modes\n')
        fprintf('is decremented by 1 \n')
        mod_cut=mod_cut+1;
    elseif length(pos)<idpar.discard_threshold_factor*(ndim+1)
        % cut the points
        fprintf('\n WARNING: Mode %d has less than %d times the number of of PV entries (i.e. %d) \n',i,idpar.discard_threshold_factor, ndim+1 )
        fprintf(' Its points are discarded and the number of modes\n')
        fprintf(' is decremented by 1 \n')
        mod_cut=mod_cut+1;
    else
        ind_inl_new=[ind_inl_new;inliers.pos(pos)];
        % note the use of mod_cut in the next line:
        % if an inlier was in the cluster 4 and 2 clusters have been cut
        % it is re-attributed to the cluster 4-2=2.
        ind_class_inl_new=[ind_class_inl_new;inliers.class(pos)-mod_cut];
    end
end
% Update the number of models
s_new=s-mod_cut;
if mod_cut>0
    fprintf('\n The new number of modes after removing %d mode(s)',mod_cut)
    fprintf('\n is: %d \n',s_new)
    if s_new<=0
        fprintf('\n ======================================================== ');
        fprintf('\n REGRESSION ERROR: according to the criteria for'); 
        fprintf('\n discarding modes, no mode survived !');
        fprintf('\n hit_regression.m ends here. ');
        fprintf('\n ======================================================== \n');
        error('hit_regression:end','hit_regression ends here.')
    end
end
inliers_new.pos=ind_inl_new;
inliers_new.class=ind_class_inl_new;