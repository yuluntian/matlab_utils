% function [R, info] = rotation_averaging_gauss_newton(measurements, R, options)
% A simple Gauss-Newton solver for 2D/3D rotation averaging.
% 
% Yulun Tian
function [R, info] = rotation_averaging_gauss_newton(measurements, R, options)
%fprintf('=== Begin Rotation Averaging Gauss-Newton ===\n\n');
if nargin < 3
    options = struct;
end
if ~isfield(options, 'tangent_space_parametrization')
    options.tangent_space_parametrization = 'global';
end
if ~isfield(options, 'rotation_distance')
    options.rotation_distance = 'chordal';
end
if ~isfield(options, 'lambda')
    options.lambda = 0;
end
if ~isfield(options, 'max_iterations')
    options.max_iterations = 100;
end
if ~isfield(options, 'gradnorm_tol')
    options.gradnorm_tol = 1e-2;
end
if ~isfield(options, 'verbose')
    options.verbose = true;
end

info = struct;
d = size(measurements.R{1},1);
n = max(max(measurements.edges));
% If using chordal distance, construct additional problem data for faster
% computation
problem_data = struct;
if strcmp(options.rotation_distance, 'chordal')
    problem_data_start = tic;
    print_if(options.verbose, 'Constructing problem data... ');
    problem_data.ConLap = construct_connection_Laplacian(measurements);
    problem_data.Manifold = stiefelstackedfactory(n, d, d);
    info.problem_data_time = toc(problem_data_start);
    print_if(options.verbose, 'Constructing problem data time: %.1f sec\n', info.problem_data_time);
end

% Save optimization stats
info.iterations = [];
info.costs = [];
info.gradnorms = [];
info.linearize_time = [];
info.solve_time = [];
if isfield(options, 'eval_func')
    info.eval_results = [];
end
iter = 0;
while true
    
    % Evaluate cost and gradients
    cost = evaluate_rotation_averaging_cost(measurements, R, options, problem_data);
    g = differentiate_rotation_averaging(measurements, R, options, problem_data);
    gradnorm = norm(g);
    info.iterations(end+1) = iter;
    info.costs(end+1) = cost;
    info.gradnorms(end+1) = gradnorm;
    if isfield(options, 'eval_func')
        info.eval_results = [info.eval_results options.eval_func(R)];
    end
    if iter >= options.max_iterations || gradnorm < options.gradnorm_tol
        break;
    end

    % Linearize and construct Gauss-Newton system
    linearize_start = tic;
    [r, J, W] = linearize_rotation_averaging(measurements, R, options);
    grad = J' * (W * r);
    H = J' * (W * J) + options.lambda * speye(size(J,2));
    info.linearize_time(end+1) = toc(linearize_start);
    
    % Solve Gauss-Newton system
    solve_start = tic;
    x = - H \ grad;
    info.solve_time(end+1) = toc(solve_start);
    
    % Apply tangent space solution
    R = rotation_averaging_exp(R, x, options);
    if options.verbose
        fprintf('Iter=%i, cost=%f, gradnorm=%.2e, xnorm=%.2e \n', ...
                 iter, cost, gradnorm, norm(x));
    end
    iter = iter + 1;
end
if options.verbose
    fprintf('Final result: iter=%i, cost=%f, gradnorm=%.2e. \n', ...
                       iter, cost, gradnorm);
end
if gradnorm > info.gradnorms(1)
    warning('Terminate with larger gradient norm: %.3e vs init %.3e', ...
        gradnorm, info.gradnorms(1));
end
%fprintf('=== End Rotation Averaging Gauss-Newton ===\n\n');
end