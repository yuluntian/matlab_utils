clear; clc; close all;

% by simulating a noiseless problem
datadir = '/home/yulun/git/dc2_pgo/data';
dataset = 'rim.g2o';
g2o_file = fullfile(datadir, dataset);
measurements = load_g2o_data(g2o_file);

% Use simulation
% sim_opts.num_rows = 5;
% sim_opts.num_cols = 5;
% sim_opts.num_floors = 5;
% sim_opts.t_stddev = 0;
% sim_opts.deg_stddev = 0;
% [measurements, true_pose, gt_info] = simulate_single_grid_pgo(sim_opts);
% Rgt = [true_pose.R{:}];

d = size(measurements.R{1}, 1);
assert(d == 2 || d == 3);
if d == 2
    p = 1;
else
    p = 3;
end
n = max(max(measurements.edges));
m = size(measurements.edges, 1);

% Use chordal initialization
[Rchordal, ~] = chordal_initialization(measurements);

% Test implementation of Riemannian Newton in matlab_utils
newton_opts = struct;
newton_opts.rotation_distance = 'chordal';
newton_opts.tangent_space_parametrization = 'local';
newton_opts.gradnorm_tol = 1e-4;
newton_opts.lambda = 0;
newton_opts.max_iterations = 50;
newton_opts.quotient_optimization = true;
[RNewton, info_newton] = rotation_averaging_newton(measurements, Rchordal, newton_opts);

% Test verification on correct solution
fprintf('\nNewton: \n')
verify_opts = struct;
verify_opts.lanczos_maxit = 1000;
verify_opts.lanczos_subspace_dim = 100;
verify_opts.verbose = true;
[is_optimal, verify_info] = certify_chordal_rotation_averaging(measurements, RNewton, verify_opts);
fprintf('is_optimal: %i, lambda_min: %.3e, max multiplier symerror: %.3e\n', is_optimal, verify_info.lambda_min, max(verify_info.multiplier_symmetric_errors));

% Test gradient 
g = differentiate_rotation_averaging(measurements, RNewton);
Mg = reshape(g, p, n);





