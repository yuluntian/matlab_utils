clear; clc; close all;

G = graph(bucky);
L = laplacian(G);
n = size(L,1);
p = 3;

nv = ones(n,1);
nv = nv / norm(nv);
for trial = 1:100
    B = randn(n,p);
    B = B - nv * (nv'*B);
    X = laplacian_solver(L,B);
    resnorm = norm(L*X-B, 'fro');
    assert(resnorm < 1e-10);
end

fprintf('Ok.\n');
