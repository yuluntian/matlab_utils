clear; clc; close all;
for trial = 1:1000
    w1 = randn(3,1);
    R1 = exp3d(w1);
    w2 = log3d(R1);
    R2 = exp3d(w2);
    err = norm(R1 - R2, 'fro');
    assert(err < 1e-6);
end
fprintf('Ok.\n');