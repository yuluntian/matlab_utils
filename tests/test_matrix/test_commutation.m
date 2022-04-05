clear; clc; close all;

m = 10;
n = 20;
K = commutation_matrix(m,n);

for trial = 1:1000
A = randn(m,n);
lhs = K * mat2vec(A);
rhs = mat2vec(A');
err = norm(lhs - rhs);
assert(err < 1e-8)
end
fprintf('Ok.\n')