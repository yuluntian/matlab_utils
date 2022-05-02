% function [R,t] = chordal_pgo_gauss_newton_2d(measurements, R, t, options)
% A simple Gauss-Newton solver for 2D pose graph optimization.
% 
% Yulun Tian
function [R, t] = pgo_gauss_newton(measurements, R, t, options)
fprintf('=== Begin PGO Gauss-Newton ===\n\n');
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
    options.lambda = 1e1;
end
if ~isfield(options, 'max_iterations')
    options.max_iterations = 100;
end
if ~isfield(options, 'gradnorm_tol')
    options.gradnorm_tol = 1e-2;
end

iter = 1;
while true
    [r, J, W] = linearize_pgo(measurements, R, t, options);
    cost = r' * (W * r);
    grad = J' * (W * r);
    gradnorm = norm(grad);
    if iter > options.max_iterations
        break;
    end
    if gradnorm < options.gradnorm_tol
        break;
    end
    % Solve Gauss-Newton system
    H = J' * (W * J) + options.lambda * speye(size(J,2));
    x = - H \ grad;
    % Apply tangent space solution
    [R, t] = pgo_exp(R, t, x, options);
    fprintf('Iter=%i, cost=%f, gradnorm=%.2e, xnorm=%.2e \n', ...
              iter, cost, gradnorm, norm(x));
    iter = iter + 1;
end

fprintf('Final result: iter=%i, cost=%f, gradnorm=%.2e. \n', ...
                   iter, cost, gradnorm);
fprintf('=== End PGO Gauss-Newton ===\n\n');
end