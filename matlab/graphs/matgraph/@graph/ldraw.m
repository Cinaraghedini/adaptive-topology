function ldraw(g,line_style)
% ldraw(g,line_style) --- draw a graph with vertices marked with their labels
% If the graph is unlabled, we use the vertex numbers instead.
% See also draw, cdraw, and ndraw.


if ~is_labeled(g)
    if nargin ==1
        ndraw(g);
        return
    else
        ndraw(g,line_style)
        return
    end
end

if nargin == 1
    draw(g);
else
    draw(g,line_style);
end

xy = getxy(g);
n = nv(g);

for v=1:n
    x = xy(v,1);
    y = xy(v,2);
    text(x,y,get_label(g,v)); 
end