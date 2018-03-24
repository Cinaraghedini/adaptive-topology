function [connected ] = graph_connected(position,param)

[A] =  initialize_matrixAdj(position,param);

g=graph;

set_matrix(g,A);

connected=isconnected(g);

free(g);

end

