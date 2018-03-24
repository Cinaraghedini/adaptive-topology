function [position connected] = graph_generator(set,options) 


g=graph;

position = initialize_position(set);

distance = squareform(pdist(position,'euclidean'));

[matrix matrixW matrixWN] = initialize_matrix(distance,set,options);

set_matrix(g,matrix);

connected = isconnected(g);

free_all;


