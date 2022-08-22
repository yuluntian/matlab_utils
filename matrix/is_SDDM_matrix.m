% function result = is_SDDM_matrix(M)
% check if the input matrix M is a SDDM matrix
% Every SSDM matrix can be uniquely written as M = L + X
% where L is a Laplacian matrix, and X is a nonnegative diagonal matrix
%
% Yulun Tian
function result = is_SDDM_matrix(M)
n = size(M, 1);
if size(M, 2) ~= n
    result = false;
    warning('Input matrix is not square.')
    return
end

% Extract non-diagonals
MDiag = spdiags(diag(M), 0, n, n);
A = MDiag - M;
if ~all(A >= 0, 'all')
    result = false;
    warning('Laplacian contains negative edge weights.')
    return;
end

% Form L+X 
degrees = sum(A, 2);
D = spdiags(degrees, 0, n, n);
L = D - A;
X = M - L;

% Check X is diagonal and nonnegative
if ~isdiag(X)
    result = false;
    warning('X = M - L is not diagonal.')
    return
end
x = diag(X);
if ~all(x >= -1e-12)
    result = false;
    warning('X = M - L does not contain nonnegative diagonal.')
    return
end

result = true;

end