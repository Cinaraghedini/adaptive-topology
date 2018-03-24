function hit_plot_clust(xi,cent,inlp,class_inl,markers,fig)
%HIT_PLOT_CLUST Plot clusters, cluster centers and outliers
%
% -------------------------------------------------------------------------
% DESCRIPTION
% -------------------------------------------------------------------------
%
% hit_plot_clust(xi,cent,inlp,class_inl,markers,fig)
% 
% -------------------------------------------------------------------------
% INPUT
% -------------------------------------------------------------------------
% xi: matrix of xi-points to be plot. Each point is a row
% 
% cent(i:): center of the i-th cluster.
%
% inlp: indexes of xi-points that are inliers. If inl(j)=i this means that
% xi(i,:) is an inlier. 
%
% class_inl: classification of the inliers. If class_inl(j)=i this means
% that xi(inl(j),:) is in the i-th cluster.
%
% markers: string of colored markers for plotting points belonging to
% different clusters.
%
% fig: (optional) # of the figure displaying the plot.
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

if nargin>=6
    figure(fig);
else 
    figure;
end
clf;

%number of clusters
nc=size(cent,1);

% allows to plot only points without the center
% useful for plotting unclassified FVs

if nc==0 &~isempty(class_inl)
    nc1=1;
    nc2=0;
else
    nc1=nc;
    nc2=nc;
end

% append markers if necessary
if (length(markers)<nc)
    fprintf('\n hit_plot_clust3d: cluster markers replicated for plotting (plotpar.marker_modes was not long enough)');
end
while(length(markers)<nc)
    markers=[markers;markers];
end
%number of FV
nFV=size(xi,1);

%compute the positions of the outliers
outlp=setdiff([1:nFV],inlp);
outl=xi(outlp);
%number of outliers
no=length(outl);

hold on


% dimension of the FVs
dim=size(xi,2);

if dim>=3
    % Plot the outlier as squares

    for j=1:no
        plot3(outl(j,1),outl(j,2),outl(j,3),'ms', 'LineWidth',2,'MarkerSize',16);
    end

    % Plot clusters
    hp=[];
    for j=1:nc1
        %extract the inliers
        tokeep=find(class_inl==j);
        tokeep1=inlp(tokeep);
        Sj=xi(tokeep1,:);
        for i=1:size(Sj,1)
            hp(i) = plot3(Sj(i,1),Sj(i,2),Sj(i,3),markers(j,:));
            set(hp(i), 'MarkerSize', 10,'LineWidth',1);
        end
    end

    % Plot centers as red crosses
    for j=1:nc2
        plot3(cent(j,1),cent(j,2),cent(j,3),'+r', 'LineWidth',2.5,'MarkerSize', 12);
    end
    grid


    axis square
    grid on
 
    xlabel('{(\xi^{LS , j})_1}','FontSize',16)
    ylabel('{(\xi^{LS , j})_2}','FontSize',16)
    zlabel('{(\xi^{LS , j})_3}','FontSize',16)

    %xlabel('{(\theta^{LS , j})_1}','FontSize',16)
    %ylabel('{(\theta^{LS , j})_2}','FontSize',16)
    %zlabel('{m_j}','FontSize',16)
    set(gca,'FontSize',13)
    view(3)
    hold off

else %the dimension is 2
    
    % Plot the outlier as squares

    for j=1:no
        plot(outl(j,1),outl(j,2),'ms', 'LineWidth',2,'MarkerSize',16);
    end

    % Plot clusters
    hp=[];
    for j=1:nc1
        %extract the inliers
        tokeep=find(class_inl==j);
        tokeep1=inlp(tokeep);
        Sj=xi(tokeep1,:);
        for i=1:size(Sj,1)
            hp(i) = plot(Sj(i,1),Sj(i,2),markers(j,:));
            set(hp(i), 'MarkerSize', 10,'LineWidth',1);
        end
    end

    % Plot centers as red crosses
    for j=1:nc2
        plot(cent(j,1),cent(j,2),'+r', 'LineWidth',2.5,'MarkerSize', 12);
    end
    grid


    axis square
    grid on
 
    xlabel('{(\xi^{LS , j})_1}','FontSize',16)
    ylabel('{(\xi^{LS , j})_2}','FontSize',16)
    
    %xlabel('{(\theta^{LS , j})_1}','FontSize',16)
    %ylabel('{(\theta^{LS , j})_2}','FontSize',16)
    set(gca,'FontSize',13)
    hold off
end
