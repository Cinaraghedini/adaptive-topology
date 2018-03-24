function g = graph(n)
% graph: constructor for the graph class
% graph(n) --- create a new graph handle for a graph with n vertices
% graph(h) --- copy h into a new graph
% graph    --- create an empty graph
% graph(edge_list) --- create graph from n-by-2 list of edges


% make sure system has been initialized

if ~graph_system_exists
   graph_init;
end

global GRAPH_MAGIC;

% get a slot to store this graph
idx = find_available;
if (idx == 0)
    error('Graph system memory is full; cannot create new graph');
end

% default number of vertices is 0
if (nargin==0)
    n = 0;
end

% if n is a graph, then we perform a copy ...

if (isa(n,'graph'))
    GRAPH_MAGIC.in_use(idx) = 1;
    GRAPH_MAGIC.graphs{idx} = GRAPH_MAGIC.graphs{n.idx};
    g.idx = idx;
    g = class(g,'graph');
    make_logical(g);
    return
end

% ... otherwise we create a new graph

[nr,nc] = size(n);
if (nc > 2)
    error('Graph constructor argument is wrong shape')
end
if (nc==2)
    v = max(max(n));  
    GRAPH_MAGIC.in_use(idx) = 1;
    GRAPH_MAGIC.graphs{idx}.array = logical(sparse([],[],[],v,v,1));
     g.idx = idx;
    g = class(g,'graph');
    
    if (v <= GRAPH_MAGIC.large_size) 
        full(g);
    end
    add(g,n);
    return
end
    


GRAPH_MAGIC.in_use(idx) = 1;
GRAPH_MAGIC.graphs{idx}.array = logical(sparse([],[],[],n,n,n));

g.idx = idx;
g = class(g,'graph');


% if the graph is "small" enough, convert to full storage

if (n < GRAPH_MAGIC.large_size) 
    full(g)
end

make_logical(g);