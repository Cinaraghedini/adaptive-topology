function hit_plotpar_check(lfields)
%HIT_PLOTPAR_CHECK Check if the cell array of fields lfields in plotpar is well defined
%
% -------------------------------------------------------------------------
% DESCRIPTION
% -------------------------------------------------------------------------
%
% hit_plotpar_check(lfields)
% 
% -------------------------------------------------------------------------
% INPUT
% -------------------------------------------------------------------------
% lfields: cell array of field names
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

global plotpar ;
if  ~isstruct(plotpar),
    hit_error;
end
nfields=length(lfields);
for i=1:nfields
    if strcmp(lfields{i},'marker_modes') & (~isfield(plotpar,'marker_modes') | ~iscellstr({plotpar.marker_modes(:)}))
        fprintf('\n ======================================================== ');
        fprintf('\n PLOTPAR CHECK ERROR: markers for points belonging to');
        fprintf('\n different clusters/modes are not avaliable');
        fprintf('\n Set correctly plotpar.marker_modes and retry.');
        fprintf('\n ======================================================== \n');
        error('hit_plotpar_check:end','Error in the variable plotpar.')
    end
    if strcmp(lfields{i},'color_regions') & ~isfield(plotpar,'color_regions')
        fprintf('\n ======================================================== ');
        fprintf('\n PLOTPAR CHECK ERROR: colors for different regions');
        fprintf('\n are not avaliable');
        fprintf('\n Set correctly plotpar.color_regions and retry.');
        fprintf('\n ======================================================== \n');
        error('hit_plotpar_check:end','Error in the variable plotpar.')
    end
    if strcmp(lfields{i},'color_surface') & ~isfield(plotpar,'color_surface')
        fprintf('\n ======================================================== ');
        fprintf('\n PLOTPAR CHECK ERROR: a color for plotting the PWA map');
        fprintf('\n has been not assigned');
        fprintf('\n Set correctly plotpar.color_surface and retry.');
        fprintf('\n ======================================================== \n');
        error('hit_plotpar_check:end','Error in the variable plotpar.')
    end
    if strcmp(lfields{i},'plot_xi_points_yn') & ~isfield(plotpar,'plot_xi_points_yn')
        fprintf('\n ======================================================== ');
        fprintf('\n PLOTPAR CHECK ERROR: should xi-points be plot ');
        fprintf('\n before clustering ?');
        fprintf('\n Set plotpar.plot_xi_points_yn to Y or N.');
        fprintf('\n ======================================================== \n');
        error('hit_plotpar_check:end','Error in the variable plotpar.')
    else switch lower(plotpar.plot_xi_points_yn)
            case{'y','n'}
            otherwise
                fprintf('\n ======================================================== ');
                fprintf('\n PLOTPAR CHECK ERROR: the value of');
                fprintf('\n plotpar.plot_xi_points_yn must be Y or N.');
                fprintf('\n ======================================================== \n');
                error('hit_plotpar_check:end','Error in the variable plotpar.')
        end
    end
    if strcmp(lfields{i},'plot_xi_points_fig') & lower(plotpar.plot_xi_points_yn)=='y'& (~isfield(plotpar,'plot_xi_points_fig') | isempty(plotpar.plot_xi_points_fig))
        fprintf('\n ======================================================== ');
        fprintf('\n PLOTPAR CHECK ERROR: no figure number has been specified');
        fprintf('\n for plotting xi-points before clustering.');
        fprintf('\n Set correctly plotpar.plot_xi_points_fig.');
        fprintf('\n ======================================================== \n');
        error('hit_plotpar_check:end','Error in the variable plotpar.')
    end
    if strcmp(lfields{i},'plot_class_xi_points_yn') & ~isfield(plotpar,'plot_class_xi_points_yn')
        fprintf('\n ======================================================== ');
        fprintf('\n PLOTPAR CHECK ERROR: should classified xi-points be plot ');
        fprintf('\n after clustering ?');
        fprintf('\n Set plotpar.plot_class_xi_points_yn to Y or N.');
        fprintf('\n ======================================================== \n');
        error('hit_plotpar_check:end','Error in the variable plotpar.')
    else switch lower(plotpar.plot_class_xi_points_yn)
            case{'y','n'}
            otherwise
                fprintf('\n ======================================================== ');
                fprintf('\n PLOTPAR CHECK ERROR: the value of');
                fprintf('\n plotpar.plot_class_xi_points_yn must be Y or N.');
                fprintf('\n ======================================================== \n');
                error('hit_plotpar_check:end','Error in the variable plotpar.')
        end
    end
    if strcmp(lfields{i},'plot_class_xi_points_fig') & lower(plotpar.plot_class_xi_points_yn)=='y'& (~isfield(plotpar,'plot_class_xi_points_fig') | isempty(plotpar.plot_class_xi_points_fig))
        fprintf('\n ======================================================== ');
        fprintf('\n PLOTPAR CHECK ERROR: no figure number has been specified');
        fprintf('\n for plotting classified xi-points after clustering.');
        fprintf('\n Set correctly plotpar.plot_class_xi_points_fig.');
        fprintf('\n ======================================================== \n');
        error('hit_plotpar_check:end','Error in the variable plotpar.')
    end
    if strcmp(lfields{i},'datapoints_yn') & ~isfield(plotpar,'datapoints_yn')
        fprintf('\n ======================================================== ');
        fprintf('\n PLOTPAR CHECK ERROR: should datapoints be plot ');
        fprintf('\n together with the modes ?');
        fprintf('\n Set plotpar.datapoints_yn to Y or N.');
        fprintf('\n ======================================================== \n');
        error('hit_plotpar_check:end','Error in the variable plotpar.')
    else switch lower(plotpar.datapoints_yn)
            case{'y','n'}
            otherwise
                fprintf('\n ======================================================== ');
                fprintf('\n PLOTPAR CHECK ERROR: the value of');
                fprintf('\n plotpar.datapoints_yn must be Y or N.');
                fprintf('\n ======================================================== \n');
                error('hit_plotpar_check:end','Error in the variable plotpar.')
        end
    end
end %end for
