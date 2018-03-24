function hit_plot_idmodes(X_or, y_or,inliers,idmodes,fig,flag_bigregionsim)
%HIT_PLOT_IDMODES Plot the identified PWA model if the domain is 1D or 2D.
%
% -------------------------------------------------------------------------
% DESCRIPTION
% -------------------------------------------------------------------------
%
% hit_plot_idmodes(X_or, y_or,inliers,idmodes,fig,flag_bigregionsim)
% 
% -------------------------------------------------------------------------
% INPUT
% -------------------------------------------------------------------------
% X_or: matrix containing the regressors. Each row is a datapoint.
%
% y_or: column vector containing the output datapoints.
%
% inliers: structure containing information on the inliers. Type 'help
% hit_regression' for a description of the fields.
%
% idmodes: structure containing information on the identified PWA model.
% Type 'help hit_regression' for a description of the fields.
%
% fig: # of the figure displaying the plot.
%
% flag_bigregionsim (optional, default=0): If =1, the modes over the
% regressor set for simulation (idpar.Regressor_set_sim) are plot.
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
% check the fields needed in plotpar
hit_plotpar_check({'marker_modes','color_regions','color_surface'})

if ~isfield(plotpar,'perm_colors') | isempty(plotpar.perm_colors)
    perm_colors=1:size(plotpar.color_regions,1); % no color permutations
end

if nargin<6
    flag_bigregionsim=0; %by default don't use Bigregion_sim
end

ind_class_inl=inliers.class;
ind_inl=inliers.pos;


ndim=size(X_or,2); % dimension of the domain
ndata=size(X_or,1);
m=size(idmodes.par{1},1);
s=length(idmodes.par); % n of submodels
theta_est={idmodes.par{:}}; % format of the coefficients for pwa_plot
if ndim==2
    % if no x_grid or y_grid have been supplied by the user, set both
    % automatically
    if ~isfield(plotpar,'x_grid') | ~isfield(plotpar,'y_grid') | isempty(plotpar.x_grid) | isempty(plotpar.y_grid)
        tempv=[];
        for i=1:s
            tempv=[tempv;extreme(idmodes.regions(i))];
        end
        minx=min(tempv(:,1));
        maxx=max(tempv(:,1));
        miny=min(tempv(:,2));
        maxy=max(tempv(:,2));
        % set 50 points by defaults
        plotpar.x_grid=minx:(maxx-minx)/30:maxx;
        plotpar.y_grid=miny:(maxy-miny)/30:maxy;
        clear tempv;
    end
    if ~flag_bigregionsim
        minz=hit_pwa_plot2d(theta_est,idmodes.regions,[],plotpar.x_grid,plotpar.y_grid,fig,plotpar.color_surface);
    else
        minz=hit_pwa_plot2d(theta_est,idmodes.regions_sim,[],plotpar.x_grid,plotpar.y_grid,fig,plotpar.color_surface);
    end
    hold on

    switch lower(plotpar.datapoints_yn)
        case{'y'}
            hp=zeros(ndata,1);

            % markers:  different colored markers for plotting mode datasets

            markers=plotpar.marker_modes;
            while(size(markers,1)<m)
                markers=[markers;markers];
                fprintf('Mode markers added for plotting the final modes (plotpar.marker_modes was not long enough)');
            end

            for i=1:ndata
                if isempty(setdiff([i],ind_inl))
                    hp(i) = plot3(X_or(i,1),X_or(i,2),y_or(i),markers(ind_class_inl(find(ind_inl==i)),:));
                else
                    %the i-th point is an outlier
                    hp(i) = plot3(X_or(i,1),X_or(i,2),y_or(i),'*r');
                end
                hold on
                set(hp(i), 'MarkerSize', 8,'LineWidth',1.5);
            end
    end
    colors_regions=plotpar.color_regions;

    if ~exist('displac')
        assi=axis;
        %displac=assi(5);
        displac=minz;
    end

    % permute the colors to have the same color in the same region
    mintemp=min([size(plotpar.color_regions,1),s]);
    temp=zeros(mintemp,size(plotpar.color_regions,2));
    temp=colors_regions(plotpar.perm_colors(1:mintemp),:);
    if ~flag_bigregionsim
        hit_plot_regions3d(idmodes.regions,displac,temp);
    else
        hit_plot_regions3d(idmodes.regions_sim,displac,temp);
    end

    %title('Estimated modes and classified datapoints');

    %xlabel('u(k)'),ylabel('y(k+1)')

    hold off
end

if ndim==1
    % if no x_grid has been supplied by the user, set it
    % automatically
    if ~isfield(plotpar,'x_grid') | ~isempty(plotpar.x_grid)
        tempv=[];
        for i=1:s
            tempv=[tempv;extreme(idmodes.regions(i))];
        end
        minx=min(tempv(:,1));
        maxx=max(tempv(:,1));
        % set 50 points by defaults
        plotpar.x_grid=minx:(maxx-minx)/50:maxx;
        clear tempv;
    end
    if ~flag_bigregionsim
        % minz is computed but not used ... regions are always plot at the
        % zero level
        minz=hit_pwa_plot1d(theta_est,idmodes.regions,[],fig);
    else
        minz=hit_pwa_plot1d(theta_est,idmodes.regions_sim,[],fig);
    end
    hold on
    switch lower(plotpar.datapoints_yn)
        case{'y'}
            hp=zeros(ndata,1);
            % markers:  different color/shapes for plotting points belonging to
            % 			  different clusters
            markers=plotpar.marker_modes;
            for i=1:ndata
                if isempty(setdiff([i],ind_inl))
                    hp(i) = plot(X_or(i,1),y_or(i),markers(ind_class_inl(find(ind_inl==i)),:));
                else
                    %the i-th point is an outlier
                    hp(i) = plot(X_or(i,1),y_or(i),'*y');
                end
                hold on
                set(hp(i), 'MarkerSize', 8,'LineWidth',1.5);
            end
    end
    hold off
end

if ndim>2
    disp('hit_plot_models does not create a plot: the regressor set has more than 2 dimensions')
end
