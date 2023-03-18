clear; clc; close all;

% This function requires the DC2-PGO library
dc2pgo_dir = '/home/yulun/git/dc2_pgo/matlab';
run(fullfile(dc2pgo_dir, 'setup.m'));

%%  Load a dataset
datadir = '/home/yulun/git/dc2_pgo/data';
dataset = 'cubicle.g2o';
g2o_file = fullfile(datadir, dataset);
measurements = load_g2o_data(g2o_file);

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
[R, ~] = chordal_initialization(measurements);

%% Test cost function

opts = struct;
opts.tangent_space_parametrization = 'global';
opts.rotation_distance = 'chordal';
problem_data = struct;
problem_data.ConLap = construct_connection_Laplacian(measurements);

cost1 = evaluate_rotation_averaging_cost(measurements, R, opts);
cost2 = evaluate_rotation_averaging_cost_chordal(measurements, R, opts, problem_data);
cost_err = abs(cost1 - cost2);
assert(cost_err < 1e-6);

%% Test gradient function

g1 = differentiate_rotation_averaging(measurements, R, opts);
g2 = differentiate_rotation_averaging_chordal(measurements, R, opts, problem_data);
grad_err = norm(g1-g2);
assert(grad_err < 1e-6);


