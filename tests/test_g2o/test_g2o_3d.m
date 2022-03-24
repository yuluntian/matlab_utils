clear; clc; close all;

%% Dependency: 
% Running this script requires the MATLAB implementation of 
% Tian et al, "Distributed Certifiably Correct Pose-Graph Optimization"
RBCD_root_dir = '/home/yulun/git/bcm2019/code/MATLAB/';
run(fullfile(RBCD_root_dir, 'setup.m'));

%% Generate synthetic data
sim_opts.num_rows = 12; % height of a single robot grid
sim_opts.num_columns = 12; % width of a single robot grid
sim_opts.num_floors = 12; % number of floors of a single grid
sim_opts.lc_prob = 0.3; % probability of intra-robot loop closure
sim_opts.fov = 1.1; % field of view, within which inter-robot loop closure is sampled 
sim_opts.velocity = 1; % velocity of each robot
sim_opts.t_stddev = 1e-1; % translational measurement standard deviation (meter)
sim_opts.deg_stddev = 1.5; % rotational measurement standard deviation (degree)
sim_opts.lc_outlier_prob = 0.0; % percentage of outlier loop closures

% simulate! 
[measurements, true_pose, gt_info] = simulate_single_grid_3D(sim_opts);
d = 3;
n = max(max(measurements.edges));

%% Obtain initial guess
[R, t] = chordal_initialization(measurements);
xhat.R = R;
xhat.t = t;

%% IO
filename = 'test.g2o';
write_to_g2o(filename, measurements, xhat);
[measurements2, poses] = load_from_g2o(filename);

%% Sanity checks
tol = 1e-5;

% Check vertices
assert(all(poses.vertices' == 1:n));  % Vertex IDs should be consecutive
R2 = [poses.R{:}];
t2 = [poses.t{:}];
for i = 1: n
    R1 = xhat.R(:, (i-1)*d+1:i*d);
    R2 = poses.R{i};
    t1 = xhat.t(:, i);
    t2 = poses.t{i};
    assert(norm(R1 - R2, 'fro') < tol);
    assert(norm(t1 - t2, 'fro') < tol);
end

% Check edges
for k = 1:size(measurements.edges,1)
    
    assert(measurements.edges(k,1) == measurements2.edges(k,1));
    assert(measurements.edges(k,2) == measurements2.edges(k,2));
    
    dR = measurements.R{k};
    dt = measurements.t{k};
    kappa = measurements.kappa{k};
    tau = measurements.tau{k};
    
    dR2 = measurements2.R{k};
    dt2 = measurements2.t{k};
    kappa2 = measurements2.kappa{k};
    tau2 = measurements.tau{k};
    
    assert(norm(dR - dR2, 'fro') < tol);
    assert(norm(dt - dt2) < tol);
    assert(abs(kappa - kappa2) < tol);
    assert(abs(tau - tau2) < tol);
    
end

fprintf('Ok.\n');




