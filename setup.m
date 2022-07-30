%% Dependency
manopt_dir = '/home/yulun/git/Manopt_7.0/manopt';
manopt_import_file = fullfile(manopt_dir, 'importmanopt.m');
run(manopt_import_file);

%% Add this directory
addpath(genpath(pwd));