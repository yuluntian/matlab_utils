% function [R, t, info] = pgo_gauss_newton(measurements, R, t, options)
% A simple Gauss-Newton solver for 2D/3D pose graph optimization.
% 
% Yulun Tian
function [R, t, info] = pgo_gauss_newton(measurements, R, t, options)

if nargin < 4
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

print_if(options.verbose, '=== Begin PGO Gauss-Newton ===\n\n');

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
    [r, J, W] = linearize_pgo(measurements, R, t, options);
    cost = evaluate_pgo_cost(measurements, R, t);
    g = differentiate_pgo(measurements, R, t);
    gradnorm = norm(g);
    info.iterations(end+1) = iter;
    info.costs(end+1) = cost;
    info.gradnorms(end+1) = gradnorm;
    if isfield(options, 'eval_func')
        info.eval_results = [info.eval_results options.eval_func(R, t)];
    end
    if iter >= options.max_iterations || gradnorm < options.gradnorm_tol
        break;
    end
    % Solve Gauss-Newton system
    grad = J' * (W * r);
    H = J' * (W * J) + options.lambda * speye(size(J,2));
    x = - H \ grad;
    % Apply tangent space solution
    [R, t] = pgo_exp(R, t, x, options);
    print_if(options.verbose, 'Iter=%i, cost=%f, gradnorm=%.2e, xnorm=%.2e \n', ...
              iter, cost, gradnorm, norm(x));
    iter = iter + 1;
end

print_if(options.verbose, 'Final result: iter=%i, cost=%f, gradnorm=%.2e. \n', ...
                   iter, cost, gradnorm);
print_if(options.verbose, '=== End PGO Gauss-Newton ===\n\n');
end