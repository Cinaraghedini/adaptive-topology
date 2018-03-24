function H=hit_holes(idmodes,Bigregion)
%HIT_HOLES Find if the union of the mode regions leave some "holes" in the regressor set
% -------------------------------------------------------------------------
% DESCRIPTION
% -------------------------------------------------------------------------
%
% H=hit_holes(idmodes,Bigregion)
% 
% -------------------------------------------------------------------------
% INPUT
% -------------------------------------------------------------------------
% idmodes: structure containing all information on the PWA function.
% Type 'help hit_regression' for a description of the fields. 
%
% Bigregion: polytope describing the regression set.
%
% -------------------------------------------------------------------------
% OUTPUT                                                                                                   
% -------------------------------------------------------------------------
%
% H: polytope array equal to set difference between Bigregion and the union
% of mode regions. If it is an empty polytope , no "hole" has been left.
%

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

H=Bigregion\idmodes.regions;