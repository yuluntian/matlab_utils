% function L = construct_weighted_laplacian(nodes, edges, weights)
% Given a list of nodes and weighted edges, return the corresponding
% weighted Laplacian matrix using the given node ordering.
%
% Yulun Tian
function L = construct_weighted_laplacian(nodes, edges, weights)

n = length(nodes);
r = [];
c = [];
v = [];
assert(length(weights) == size(edges, 1));

for k = 1:size(edges,1)
    v1 = edges(k,1);
    v2 = edges(k,2);
    w = weights(k);
    i1 = find(nodes == v1);
    i2 = find(nodes == v2);
    assert(length(i1) == 1)
    assert(length(i2) == 1)
    
    r(end+1) = i1;
    c(end+1) = i1;
    v(end+1) = w;
    
    r(end+1) = i1;
    c(end+1) = i2;
    v(end+1) = -w;

    r(end+1) = i2;
    c(end+1) = i1;
    v(end+1) = -w;

    r(end+1) = i2;
    c(end+1) = i2;
    v(end+1) = w;

end

L = sparse(r, c, v, n, n);

end