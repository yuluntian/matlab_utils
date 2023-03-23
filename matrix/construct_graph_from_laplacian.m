% function [G, nodes, edges, weights] = construct_graph_from_laplacian(L, options)
% Create a matlab weighted graph object corresponding to the input weighted
% Laplacian L.
%
% Yulun Tian
function [G, nodes, edges, weights] = construct_graph_from_laplacian(L, options)
if nargin < 2
    options = struct;
end
if ~isfield(options, 'weighted')
    options.weighted = true;
end
n = size(L,1);
assert(size(L,2) == n);

% Extract edges by getting nonzeros from the upper triangular portion of L
[i, j, v] = find(triu(L, 1));
m = length(i);
if options.weighted
    weights = -v;
    G = graph(i, j, weights);
else
    weights = ones(m,1);
    G = graph(i, j);
end

nodes = 1:n;
edges = [i j];

end