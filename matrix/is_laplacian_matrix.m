% function result = is_laplacian_matrix(L)
% Return if the given square matrix L is a Laplacian
%
% Yulun Tian
function result = is_laplacian_matrix(L)

n = size(L, 1);
if size(L, 2) ~= n
    result = false;
    warning('Input matrix is not square.')
    return
end
% L = D - A
degrees = diag(L);
D = spdiags(degrees, 0, n, n);
A = D - L;

% All weights in adjacency matrix should be nonnegative
if ~all(A >= 0, 'all')
    result = false;
    warning('Laplacian contains negative edge weights.')
    return;
end

% Diagonal should contain sum of edge weights for each vertex
tmp = sum(L, 2);
if norm(tmp) > 1e-10
    result = false;
    warning('Row sum does not equal zero')
    return;
end

result = true;

end