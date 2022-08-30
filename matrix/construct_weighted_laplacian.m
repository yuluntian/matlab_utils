% function L = construct_weighted_laplacian(nodes, edges, weights)
% Given a list of nodes and weighted edges, return the corresponding
% weighted Laplacian matrix using the given node ordering.
%
% Yulun Tian
function L = construct_weighted_laplacian(nodes, edges, weights)

n = length(nodes);
m = size(edges, 1);
r = zeros(1, 4*m);  % Each edge adds 4 values to the Laplacian matrix
c = zeros(1, 4*m);
v = zeros(1, 4*m);
assert(length(weights) == size(edges, 1));

for k = 1:size(edges,1)
    offset = (k-1) * 4;
    v1 = edges(k,1);
    v2 = edges(k,2);
    w = weights(k);
    i1 = find(nodes == v1);
    i2 = find(nodes == v2);
    assert(length(i1) == 1)
    assert(length(i2) == 1)
    
    r(offset+1) = i1;
    c(offset+1) = i1;
    v(offset+1) = w;
    
    r(offset+2) = i1;
    c(offset+2) = i2;
    v(offset+2) = -w;

    r(offset+3) = i2;
    c(offset+3) = i1;
    v(offset+3) = -w;

    r(offset+4) = i2;
    c(offset+4) = i2;
    v(offset+4) = w;

end

L = sparse(r, c, v, n, n);

end