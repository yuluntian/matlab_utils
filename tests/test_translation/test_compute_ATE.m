clear; clc; close all;
d = 3;
n = 1000;

for trial = 1:100
    t_src = 100 * rand(d,n);
    
    transform = struct;
    transform.R = randrot(d);
    transform.t = rand(d,1);
    
    t_dst = transform_translations(t_src, transform);
    
    ATE = compute_ATE(t_src, t_dst);
    assert(ATE < 1e-8);
    
    t_dst = t_dst + 0.01 * randn(d,n);
    ATE = compute_ATE(t_src, t_dst);
end

fprintf('Ok.\n')