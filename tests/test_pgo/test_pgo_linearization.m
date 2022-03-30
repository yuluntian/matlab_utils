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

gn_options = struct;
gn_options.rotation_distance = 'chordal';
gn_options.tangent_space_parametrization = 'global';
[rg, Jg, Wg] = linearize_pgo(measurements, R, t, gn_options);
gn_options.tangent_space_parametrization = 'local';
[rl, Jl, Wl] = linearize_pgo(measurements, R, t, gn_options);

tol = 1e-5;
assert(norm(rg - rl) < tol);
assert(norm(Wg - Wl, 'fro') < tol);