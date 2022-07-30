clear; clc; close all;
%% NOTE: this script requires the DC2-PGO repo
DC2PGO_root = '/home/yulun/git/dc2_pgo';
dataset_dir = fullfile(DC2PGO_root, 'data');
setup_file = fullfile(DC2PGO_root, 'matlab/setup.m');
run(setup_file);

%% Test vee and hat operations
SO = rotationsfactory(3);
R = SO.rand();
H = randn(3,3);
S = SO.proj(R,H);
S2 = hat3(vee3(S));
assert(norm(S-S2) < 1e-8);

%% Test functions to convert rotations between 2D matrix and 3D tensor
M = rotationsfactory(3, 5);
R = M.rand();
Rflat = rotations_tensor_to_flat(R);
R2 = rotations_flat_to_tensor(Rflat);
assert(sum(abs(R-R2), 'all') < 1e-8);

%% Test implementation of cost, (euclidean) grad and hess operators
% by simulating a noiseless problem
% g2o_file = '../data/CSAIL.g2o';
% measurements = load_g2o_data(g2o_file);

% Use simulation
sim_opts.num_rows = 5;
sim_opts.num_cols = 5;
sim_opts.num_floors = 5;
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
newton_opts.rotation_distance = 'chordal';
newton_opts.tangent_space_parametrization = 'local';
newton_opts.gradnorm_tol = 1e-6;
newton_opts.lambda = 0;
newton_opts.max_iterations = 50;
RNewton = rotation_averaging_newton(measurements, Rinit, newton_opts);
error_newton = compute_rotation_RMSE(RNewton, Rgt);

%% Test verification on correct solution
verify_opts = struct;
verify_opts.verbose = false;
[is_optimal, verify_info] = certify_chordal_rotation_averaging(measurements, RNewton, verify_opts);
fprintf('\nNewton: \n')
fprintf('is_optimal: %i, lambda_min: %.3e, max multiplier symerror: %.3e\n', is_optimal, verify_info.lambda_min, max(verify_info.multiplier_symmetric_errors));

%% Test verification on suboptimal solutions
% Use random initialization
M = rotationsfactory(3, size(R0,3));
Rrand = rotations_tensor_to_flat(M.rand());
[is_optimal, verify_info] = certify_chordal_rotation_averaging(measurements, Rrand, verify_opts);
fprintf('\nRandom: \n')
fprintf('is_optimal: %i, lambda_min: %.3e, max multiplier symerror: %.3e\n', is_optimal, verify_info.lambda_min, max(verify_info.multiplier_symmetric_errors));

fprintf('\nOk.\n')










