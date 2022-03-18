clear; clc; close all;

for trial = 1:100
    % Simulate a random set of rotations
    d = 3;
    n = 10;
    R_src = [];
    for i = 1:n
        R_src = [R_src quat2rotm(randrot)];
    end
    
    % Simulate a random global rotation
    R_dst_src = quat2rotm(randrot);
    R_dst = R_dst_src * R_src;
    
    [R_dst_src_est, R_src_aligned] = align_rotations(R_src, R_dst);
    error = norm(R_dst_src_est - R_dst_src, 'fro');
    assert(error < 1e-5);
end

fprintf('Ok.\n');