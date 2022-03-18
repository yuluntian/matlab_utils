clear; clc; close all;

for trial = 1:1000
    x = randn(6,1);
    X = pose_hat3(x);
    T = expm(X);
    T2 = pose_exp3(x);
    error = norm(T - T2, 'fro');
    assert(error < 1e-5);
end

fprintf('Ok.\n');