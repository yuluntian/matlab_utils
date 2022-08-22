clear; clc; close all;
%% NOTE: this script requires the DC2-PGO repo
DC2PGO_root = '/home/yulun/git/dc2_pgo';
dataset_dir = fullfile(DC2PGO_root, 'data');
setup_file = fullfile(DC2PGO_root, 'matlab/setup.m');
run(setup_file);

%% Test implementation of cost, (euclidean) grad and hess operators
% by simulating a noiseless problem
% g2o_file = '../data/CSAIL.g2o';
% measurements = load_g2o_data(g2o_file);

% Use simulation
sim_opts.num_rows = 10;
sim_opts.num_cols = 10;
sim_opts.num_floors = 10;
sim_opts.t_stddev = 0;
sim_opts.deg_stddev = 0;
[measurements, true_pose, gt_info] = simulate_single_grid_pgo(sim_opts);
Rgt = [true_pose.R{:}];

% Use perturbed chordal initialization
[Rchordal, ~] = chordal_initialization(measurements);
R0 = rotations_flat_to_tensor(Rchordal);
% Perturb initial rotation by small noise
for i = 1:size(R0,3)
    dR = Langevin_sampler_3D(eye(3), 1e2);
    R0(:,:,i) = R0(:,:,i) * dR;
end
Rinit = rotations_tensor_to_flat(R0);
error_init = compute_rotation_RMSE(Rinit, Rgt);

%% Test implementation of Riemannian Newton in matlab_utils
newton_opts = struct;
newton_opts.rotation_distance = 'geodesic';
newton_opts.tangent_space_parametrization = 'global';
newton_opts.gradnorm_tol = 1e-6;
newton_opts.lambda = 0;
newton_opts.max_iterations = 50;
newton_opts.quotient_optimization = true;
RNewton = rotation_averaging_newton(measurements, Rinit, newton_opts);
error_newton = compute_rotation_RMSE(RNewton, Rgt);

fprintf('\nOk.\n')










