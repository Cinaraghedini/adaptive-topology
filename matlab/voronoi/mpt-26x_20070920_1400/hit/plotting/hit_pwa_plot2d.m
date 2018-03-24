function minz=hit_pwa_plot2d(theta,regions,regressor_set,x1,x2,fig,color_surface)
%HIT_PWA_PLOT2D Plot a PWA function with 2D domain.
%
% -------------------------------------------------------------------------
% DESCRIPTION
% -------------------------------------------------------------------------
%
% minz=hit_pwa_plot2d(theta,regions,regressor_set,x1,x2,fig,color_surface)
% 
% -------------------------------------------------------------------------
% INPUT
% -------------------------------------------------------------------------
% theta{i}: PV of the i-th mode.
%
% regions(i): polyhedron describing the region of the i-th mode. 
%
% regressor_set: polytope specifying the polytopic domain of the function.
% It can be empty if the regions are already bounded.
%
% nfig: # of the figure displaying the plot.
%
% x_lab, y_lab: labels for the x-y axis.
%
% x1,x2: grids of points in the x,y directions, respectively.
%
% fig: number of the figure in which the plot appears.
%
% color_surface: RGB vector defining the color of the PWA surface.
% 
% -------------------------------------------------------------------------
% OUTPUT                                                                                                   
% -------------------------------------------------------------------------
% minz: minimum of the PWA function.

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

nreg=length(theta); % # of regions and models

[X1,X2] = meshgrid(x1,x2);

lx1=length(x1);
lx2=length(x2);

Z=zeros(lx2,lx1);

for i=1:lx1
    for j=1:lx2
        point=[x1(i),x2(j)]';
        Z(j,i)=hit_pwa(theta,regions,point);
    end
end
minz=min(min(Z));
figure(fig)
clf
colormap(color_surface)
mesh(X1,X2,Z,100*ones(size(X1)),'FaceColor','none') %0*ones(size(X1))
grid on
hold on
% ylabel('{u(k-1)}','FontSize',16)
ylabel('{x_2}','FontSize',16)
% xlabel('{y(k-1)}','FontSize',16)
xlabel('{x_1}','FontSize',16)
% zlabel('{y(k)}','FontSize',16)
zlabel('{y}','FontSize',16)
set(gca,'FontSize',13)
axis square
