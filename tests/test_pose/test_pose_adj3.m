clear; clc; close all;

for trial = 1:1000
    x = randn(6,1);
    T = pose_exp3(x);
    Ad = pose_adj3(T);
    
    u = randn(6,1);
    lhs =  expm(pose_hat3(Ad * u));
    rhs = T * pose_exp3(u) * inv(T);
    
    error = norm(lhs-rhs, 'fro');
    assert(error < 1e-8);
end

fprintf('Ok.\n')