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

% Should be symmetric
sym_err = norm(L - L', 'fro') / n;
if sym_err > 1e-8
    result = false;
    warning('Laplacian is not sufficiently symmetric.')
    return;
end

% All weights in adjacency matrix should be nonnegative
if ~all(A >= 0, 'all')
    result = false;
    warning('Laplacian contains negative edge weights.')
    return;
end

% Diagonal should contain sum of edge weights for each vertex
tmp = sum(L, 2);
row_sum_err = norm(tmp);
if row_sum_err/n > 1e-8
    result = false;
    warning('Row sum does not equal zero: %.3e', row_sum_err)
    return;
end

result = true;

end