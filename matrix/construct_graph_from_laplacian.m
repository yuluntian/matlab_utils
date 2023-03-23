% function [G, nodes, edges, weights] = construct_graph_from_laplacian(L)
% Create a matlab weighted graph object corresponding to the input weighted
% Laplacian L.
%
% Yulun Tian
function [G, nodes, edges, weights] = construct_graph_from_laplacian(L)
n = size(L,1);
assert(size(L,2) == n);

% Extract edges by getting nonzeros from the upper triangular portion of L
[i, j, v] = find(triu(L, 1));
weights = -v;

G = graph(i, j, weights);
nodes = 1:n;
edges = [i j];

end