clear; clc; 

%% This script requires using SE-Sync
% Modify below to DC2-PGO directory
dc2pgo_path = '~/git/dc2_pgo/matlab';
addpath(genpath(dc2pgo_path));

%% Main
dataset_dir = '/home/yulun/git/dc2_pgo/data';
g2o_file = fullfile(dataset_dir, 'input_M3500_g2o.g2o');
measurements = load_g2o_data(g2o_file);
[R,t] = chordal_initialization(measurements);

gn_options = struct;
gn_options.tangent_space_parametrization = 'global';
gn_options.lambda = 1e-4;
gn_options.gradnorm_tol = 1e-4;
chordal_pgo_gauss_newton_2d(measurements, R, t, gn_options);