% function X = laplacian_solver(L, B)
% Find a solution (out of many) to a Laplacian system in matrix form LX = B
% This implementation assumes that the underlying graph is connected
%
% Yulun Tian
function X = laplacian_solver(L, B)
n = size(L,1);
p = size(B,2);
assert(size(L,2) == n);
assert(size(B,1) == n);

G = construct_graph_from_laplacian(L);
[~, binsizes] = conncomp(G);
assert(length(binsizes) == 1, 'Input graph is not connected!');

% Solve by fixing first dimension of solution to zero
X = zeros(n,p);
Bc = B(2:end, :);
Lcc = L(2:end,2:end);
X(2:end,:) = Lcc \ Bc;

% Sanity check
resnorm = norm(L*X - B, 'fro');
assert(resnorm < 1e-6);

end