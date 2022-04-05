clear; clc; 

%% This script requires using SE-Sync
% Modify below to DC2-PGO directory
dc2pgo_path = '~/git/dc2_pgo/matlab';
addpath(genpath(dc2pgo_path));

%% Main
dataset_dir = '/home/yulun/git/dc2_pgo/data';
g2o_file = fullfile(dataset_dir, 'parking-garage.g2o');
measurements = load_g2o_data(g2o_file);
[R,t] = chordal_initialization(measurements);

options = struct;
options.rotation_distance = 'chordal';
options.gradnorm_tol = 1e-4;
pgo_newton(measurements, R, t, options);