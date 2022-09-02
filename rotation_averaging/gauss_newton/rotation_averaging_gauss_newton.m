% function [R, info] = rotation_averaging_gauss_newton(measurements, R, options)
% A simple Gauss-Newton solver for 2D/3D rotation averaging.
% 
% Yulun Tian
function [R, info] = rotation_averaging_gauss_newton(measurements, R, options)
fprintf('=== Begin Rotation Averaging Gauss-Newton ===\n\n');
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

% Save optimization stats
info = struct;
info.iterations = [];
info.costs = [];
info.gradnorms = [];
if isfield(options, 'eval_func')
    info.eval_results = [];
end
iter = 0;
while true
    [r, J, W] = linearize_rotation_averaging(measurements, R, options);
    cost = evaluate_rotation_averaging_cost(measurements, R, options);
    g = differentiate_rotation_averaging(measurements, R, options);
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
    % Solve Gauss-Newton system
    grad = J' * (W * r);
    H = J' * (W * J) + options.lambda * speye(size(J,2));
    x = - H \ grad;
    % Apply tangent space solution
    R = rotation_averaging_exp(R, x, options);
    fprintf('Iter=%i, cost=%f, gradnorm=%.2e, xnorm=%.2e \n', ...
             iter, cost, gradnorm, norm(x));
    iter = iter + 1;
end

fprintf('Final result: iter=%i, cost=%f, gradnorm=%.2e. \n', ...
                   iter, cost, gradnorm);
fprintf('=== End Rotation Averaging Gauss-Newton ===\n\n');
end