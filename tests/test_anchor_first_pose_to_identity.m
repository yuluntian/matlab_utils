clear; clc; close all;

num_poses = 25;

for trial = 1:1000
    % randomly generate a set of poses
    R0 = [];
    t0 = [];
    for i = 1:num_poses
        R0 = [R0 randrot(3)];
        t0 = [t0 randn(3,1)];
    end
    
    [R, t] = anchor_first_pose_to_identity(R0, t0);
    
    % checks
    assert(norm(R(:, 1:3) - eye(3), 'fro') < 1e-8);
    assert(norm(t(:,1)) < 1e-8);
    
    t_dist = compute_ATE(t, t0);
    assert(t_dist < 1e-8, 'translation distance: %.2e', t_dist);
    
    deg_dist = compute_rotation_RMSE(R, R0);
    assert(deg_dist < 1e-4, 'rotation distance: %.2e deg', deg_dist);
    
end

fprintf('Ok.\n')