function cdraw(g,coloring,line_style)
% cdraw(g,coloring) -- draw g with a given vertex coloring
% If no coloring is specified, the default is 'greedy'. 
% cdraw(g,coloring,line_style) --- lines have given line_style
% See also draw, ldraw, and ndraw.
%
% Author: Brian Towne

% edit these to change the colors 
edge_color = 'k';
vertex_color = 'k';
vertex_fill = ['k','w','b','g','r','c','m','y']; % only supports 8 colors
r = 0.15;

if nargin < 3
    line_style = '-';
end

if nargin < 2
    coloring = color(g,'greedy');
end

n = nv(g);
n2 = nv(coloring);

if ~(n==n2)
    error('Graph and coloring must have equal number of vertices.')
    return
end

if ~hasxy(g)
    embed(g);
end

xy = getxy(g);

% First draw the edges  
for u=1:n-1
    for v=u+1:n
        if has(g,u,v)
            x = xy([u,v],1);
            y = xy([u,v],2);
            line(x,y,'Color', edge_color,'LineStyle',line_style);
        end
    end
end

% Now draw the vertices
color_classes = parts(coloring);
num_colors = size(color_classes,2);

if num_colors > size(vertex_fill,2)
    error('cdraw currently does not support more than 8 colors.')
    return
end % Need to change this

for i=1:num_colors
    color_class_size = size(color_classes{i},2);
    for j=1:color_class_size
        v = color_classes{i}(j);
        x = xy(v,1);
        y = xy(v,2);
        rectangle('Position', [x-r/2, y-r/2, r, r],...
                  'Curvature', [1 1], ...
                  'EdgeColor', vertex_color, ...
                  'FaceColor', vertex_fill(i));    
    end
end    

axis equal
axis off
