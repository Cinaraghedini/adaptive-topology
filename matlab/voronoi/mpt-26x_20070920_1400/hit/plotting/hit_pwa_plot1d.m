function minz=hit_pwa_plot1d(theta,regions,domain,nfig,x_lab,y_lab)
%HIT_PWA_PLOT1D Plot a PWA function with 1D domain.
%
% -------------------------------------------------------------------------
% DESCRIPTION
% -------------------------------------------------------------------------
%
% minz=hit_pwa_plot1d(theta,regions,domain,nfig,x_lab,y_lab)
% 
% -------------------------------------------------------------------------
% INPUT
% -------------------------------------------------------------------------
% theta{i}: PV of the i-th mode.
%
% regions(i): polyhedron describing the region of the i-th mode. 
%
% domain: polytope specifying the polytopic domain of the function. It can
% be empty if the regions are already bounded.
%
% nfig: # of the figure displaying the plot.
%
% x_lab, y_lab: labels for the x-y axis.
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

if nargin<=6
    % default labels
    y_lab='y';
    if nargin<=5
        x_lab='x';
    end
end

nreg=length(theta); % # of regions and models

%intersect the regions with the domain, if it is not empty

% FOR THE MOMENT I'm ASSUMING THAT A MATRIX (usually = []) IS AN EMPTY POLYTOPE !
% TO BE IMPROVED !!!!

if ~isnumeric(domain)
    for i=1:nreg
        regions(i)=regions(i)&domain;
    end
end

figure(nfig);clf;
hold on
for i=1:nreg
    [V,R]=extreme(regions(i));    % compute extreme points and rays of polytope P
    if size(R,1)>0
        fprintf('\n ======================================================== ');
        fprintf('\n ERROR in hit_pwa_plot1d.m: a Polytope is unbounded.')
        fprintf('\n Pass a bounded function domain to hit_pwa_plot1d');
        fprintf('\n and the problem will disappear.');
        fprintf('\n ======================================================== \n')
        error('hit_pwa_plot1d: Polytope is unbounded!'); % existence of rays means polytope is unbounded
    end
    V=sort(V);
    % Note: don't use hit_pwa for computing the values of the function
    % at regions bounday, because a boundary point belongs to more than
    % one (closed) region and hit_pwa will make an "arbitrary" choice
    % of which mode is active (i.e. according to the order modes are listed
    y1(i)=theta{i}(1)*V(1)+theta{i}(2);
    y2(i)=theta{i}(1)*V(2)+theta{i}(2);
    plot(V',[y1(i),y2(i)],'LineWidth',1)
end
% compute the minimum of the function
minz=min([y1',y2']);
hold on
plot(regions)
for ii=1:nreg
    [xc,rc]=chebyball(regions(ii));
    text(xc(1),0,int2str(ii));
end
% reset title and labels set by the plot routine of the mpt toolbox
title('');
axis square
xlabel(x_lab,'FontSize',16)
ylabel(y_lab,'FontSize',16)
set(gca,'FontSize',13)
baseline=axis;


