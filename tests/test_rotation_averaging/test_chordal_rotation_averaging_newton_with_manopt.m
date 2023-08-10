clear; clc; close all;

% This function requires the DC2-PGO library
dc2pgo_dir = '/home/yulun/git/dc2_pgo/matlab';
run(fullfile(dc2pgo_dir, 'setup.m'));

%%  Load a dataset
% datadir = '/home/yulun/git/dc2_pgo/data';
% dataset = 'smallGrid3D.g2o';
% g2o_file = fullfile(datadir, dataset);
% measurements = load_g2o_data(g2o_file);

for trial = 1:50

% Use simulation
sim_opts.num_rows = 5;
sim_opts.num_cols = 5;
sim_opts.num_floors = 5;
sim_opts.t_stddev = 0.1;
sim_opts.deg_stddev = 3;
[measurements, true_pose, gt_info] = simulate_single_grid_pgo(sim_opts);
% Rgt = [true_pose.R{:}];

d = size(measurements.R{1}, 1);
assert(d == 3);
p = d*(d-1)/2;
n = max(max(measurements.edges));
m = size(measurements.edges, 1);

% Use chordal initialization
[R, ~] = chordal_initialization(measurements);

%% Compute Newton step

opts = struct;
opts.tangent_space_parametrization = 'global';
opts.rotation_distance = 'chordal';

[g, H] = differentiate_rotation_averaging(measurements, R, opts);
v = H \ -g;

%% Check that v is the solution to the Newton system using Manopt
% Construct tangent vector in extrinsic coordinate that corresponds to v
% Reference: appendix II in sparsification paper, equation (117)-(118)
eta = zeros(d, d*n);
for i = 1:n
    idxs = (i-1)*d+1 : i*d;
    vi = v(idxs);
    Ri = R(:, idxs);
    eta(:, idxs) = Ri * hat3(Ri' * vi);
end

ConLap = construct_connection_Laplacian(measurements);
M = stiefelstackedfactory(n, d, d);
euc_grad_T = ConLap * R';
rie_grad_T = M.proj(R', euc_grad_T);  % Riemannian gradient (transposed) in extrinsic coordinate
G = rie_grad_T';

euc_Hv_T = ConLap * eta';
rie_Hv_T = M.ehess2rhess(R', euc_grad_T, euc_Hv_T, eta');
HV = rie_Hv_T';

% If correct, we should have HV = -G
error = norm(HV+G, 'fro');
fprintf('Trial %i error: %.2e.\n', trial, error);
assert(error < 1e-8);
end

fprintf('Ok.\n')
