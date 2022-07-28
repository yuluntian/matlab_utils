% function R = rotation_averaging_newton(measurements, R, options)
% A simple implementation of the Newton method for 3D rotation averaging
function R = rotation_averaging_newton(measurements, R, options)
if nargin < 3
    options = struct;
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
if ~isfield(options, 'tangent_space_parametrization')
    options.tangent_space_parametrization = 'local';
end
d = size(measurements.R{1}, 1);
n = max(max(measurements.edges));
assert(d == 2 || d == 3, 'Rotation averaging problem must be 2D or 3D.');
if d == 2
    p = 1;
else
    p = 3;
end

for iter = 1 : options.max_iterations 
    cost = evaluate_rotation_averaging_cost(measurements, R, options);
    [grad, Hess] = differentiate_rotation_averaging(measurements, R, options);
    gradnorm = norm(grad);
    if gradnorm < options.gradnorm_tol
        fprintf('Final result: iter=%i, cost=%f, gradnorm=%.2e. \n', ...
                   iter, cost, gradnorm);
        break;
    end
    % Solve Newton system
    H = Hess + options.lambda * speye(p * n);
    x = - H \ grad;
    % Apply tangent space solution
    R = rotation_averaging_exp(R, x, options);
    fprintf('Iter=%i, cost=%f, gradnorm=%.2e, xnorm=%.2e \n', ...
              iter, cost, gradnorm, norm(x));
end

end