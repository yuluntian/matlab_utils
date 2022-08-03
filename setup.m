%% Dependency
manopt_dir = '/home/yulun/git/Manopt_7.0/manopt';
manopt_import_file = fullfile(manopt_dir, 'importmanopt.m');
if isfile(manopt_import_file)
    run(manopt_import_file);
else
    fprintf('Manopt not found. Some functions might not work. \n');
end

%% Add this directory
addpath(genpath(pwd));