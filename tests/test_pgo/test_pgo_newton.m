clear; clc; 

%% This script requires using SE-Sync
% Modify below to DC2-PGO directory
dc2pgo_path = '~/git/dc2_pgo/matlab';
addpath(genpath(dc2pgo_path));

%% Main
dataset_dir = '/home/yulun/git/dc2_pgo/data';
g2o_file = fullfile(dataset_dir, 'smallGrid3D.g2o');
measurements = load_g2o_data(g2o_file);
[R,t] = chordal_initialization(measurements);

options = struct;
options.rotation_distance = 'chordal';
options.gradnorm_tol = 1e-4;
[Rnewton, tnewton, info_newton] = pgo_newton(measurements, R, t, options);

% SE-Sync reference
[SDPval, Yopt, xhat, Fxhat, SE_Sync_info, problem_data] = SE_Sync(measurements);

error_rot = compute_rotation_RMSE(xhat.R, Rnewton);
error_trans = compute_ATE(xhat.t, tnewton);

assert(error_rot < 1e-4);
assert(error_trans < 1e-4);