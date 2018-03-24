% filename: subgrafosSize.m
% Purpose:  computes the subgraph of g and their sizes 
% Input: 
% - g - graph
% Output: 
% - subGrafo - array with the size of each graph partition

function subGrafo = subgrafosSize(g)

% parameters set
   
options = set_matlab_bgl_default();
options.istrans = 0;
options.nocheck=0;
options.unweighted=1;
options.full2sparse=1;

subgrafo = components(matrix(g),options);
nrElemSubGrafo = zeros(numel(subgrafo),1);
for i=1:numel(subgrafo)
    nrElemSubGrafo(subgrafo(i),1) = nrElemSubGrafo(subgrafo(i),1)+1;
end
idx = 0;
for i = 1:length(nrElemSubGrafo)
    if (nrElemSubGrafo(i,1) > 0)
        idx = idx + 1;
        subGrafo(idx,1) = nrElemSubGrafo(i,1);
    end
end