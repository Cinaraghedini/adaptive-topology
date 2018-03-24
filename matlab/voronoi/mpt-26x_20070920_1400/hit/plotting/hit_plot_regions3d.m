function handle=hit_plot_regions3d(Regions,offset,colors)
%HIT_PLOT_REGIONS3D Plot the regions of a PWA model with 2D domain.
%
% -------------------------------------------------------------------------
% DESCRIPTION
% -------------------------------------------------------------------------
%
% handle=hit_plot_regions3d(Regions,offset,colors)
% 
% -------------------------------------------------------------------------
% INPUT
% -------------------------------------------------------------------------
% Regions: polytope array. Each polytope MUST be bounded.
% offset: z-level defining the plane over which regions are plotted
% colors: color(i,:) is an RGB triplette
% 
% -------------------------------------------------------------------------
% OUTPUT                                                                                                   
% -------------------------------------------------------------------------
% handles: list of graphic handles for each polytope
%
% -------------------------------------------------------------------------
% ACKNOWLEDGMENTS                                                                                               
% -------------------------------------------------------------------------
% This function has been largely inspired by mpt_plotPWA.m. The mptOptions
% are used.

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

global mptOptions

% trick: define Options=[] and let the if..end set it to the default
% values.
% This is done in view of including a vector of Options
% among the function arguments ...

Options=[];
if ~isfield(Options,'abs_tol')
    Options.abs_tol=mptOptions.abs_tol;  % absolute tolerance
end
if ~isfield(Options, 'axis')        % axis for plot
    Options.axis='auto';
end
if ~isfield(Options, 'lpsolver')        % axis for plot
    Options.lpsolver=mptOptions.lpsolver;
end
if ~isfield(Options, 'extreme_solver')        % axis for plot
    Options.extreme_solver=mptOptions.extreme_solver;
end
if ~isfield(Options,'newfigure')
    Options.newfigure=mptOptions.newfigure;
end


if Options.newfigure,
   figure;  % open new figure window
else
   newplot; % get current figure (or create new figure)
end

handle=[];
maxlen=length(Regions);
    
for ii=1:maxlen,
            P=Regions(ii);
            nb=nconstr(P);
            dimP=dimension(P);
            [xc,rc]=chebyball(P);
            if rc<=Options.abs_tol
                disp('hit_plot_regions3d: Empty polytope detected!');
                continue
            end
            axis(Options.axis);

            [V,R,PA(ii)]=extreme(P,Options);
            if size(R,1)>0
                error('hit_plot_regions3d: Polytope is unbounded!');
            end

            % sort vertices in a cyclic way;
            x1=V(:,1);
            x2=V(:,2);

            ang=angle([(x1-xc(1))+(x2-xc(2))*sqrt(-1)]);
            [val,ind]=sort(ang);
            x1=x1(ind);
            x2=x2(ind);
            h=patch(x1,x2,offset*ones(size(x1)),colors(ii,:));
            handle=[handle;h];
            text(xc(1),xc(2),offset,int2str(ii));
end
