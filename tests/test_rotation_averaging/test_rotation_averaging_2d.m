clear; clc; close all;
%% NOTE: this script requires the DC2-PGO repo
DC2PGO_root = '/home/yulun/git/dc2_pgo';
dataset_dir = fullfile(DC2PGO_root, 'data');
setup_file = fullfile(DC2PGO_root, 'matlab/setup.m');
run(setup_file);

%% Load a dataset
g2o_file = fullfile(dataset_dir, 'input_M3500_g2o.g2o');
measurements = load_g2o_data(g2o_file);
n = max(max(measurements.edges));
% Use chordal initialization
[Rchordal, ~] = chordal_initialization(measurements);

%% Test implementation of Riemannian Newton in matlab_utils
newton_opts = struct;
newton_opts.rotation_distance = 'chordal';
newton_opts.tangent_space_parametrization = 'local';
newton_opts.gradnorm_tol = 1e-6;
newton_opts.lambda = 0;
newton_opts.max_iterations = 50;
newton_opts.quotient_optimization = true;
RNewton = rotation_averaging_newton(measurements, Rchordal, newton_opts);

%% Test verification on correct solution
fprintf('\nNewton: \n')
verify_opts = struct;
verify_opts.verbose = false;
[is_optimal, verify_info] = certify_chordal_rotation_averaging(measurements, RNewton, verify_opts);
fprintf('is_optimal: %i, lambda_min: %.3e, max multiplier symerror: %.3e\n', is_optimal, verify_info.lambda_min, max(verify_info.multiplier_symmetric_errors));

%% Test verification on suboptimal solutions
% Use random initialization
fprintf('\nRandom: \n')
M = rotationsfactory(2, n);
Rrand = rotations_tensor_to_flat(M.rand());
[is_optimal, verify_info] = certify_chordal_rotation_averaging(measurements, Rrand, verify_opts);
fprintf('is_optimal: %i, lambda_min: %.3e, max multiplier symerror: %.3e\n', is_optimal, verify_info.lambda_min, max(verify_info.multiplier_symmetric_errors));

fprintf('\nOk.\n')










