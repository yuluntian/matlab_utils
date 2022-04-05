% function [R,t] = pgo_newton(measurements, R, t, options)
% A simple implementation of the Newton method for 3D pose graph
% optimization.
%
% Yulun Tian
function [R,t] = pgo_newton(measurements, R, t, options)
if nargin < 4
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
d = length(measurements.t{1});
n = max(max(measurements.edges));

% Currently, this implementation only supports 3D and local tangent space
% parametrization
assert(d == 3, 'PGO problem is not 3D.');
p = 6;
options.tangent_space_parametrization = 'local';

for iter = 1 : options.max_iterations 
    cost = evaluate_pgo_cost(measurements, R, t, options);
    [grad, Hess] = differentiate_pgo(measurements, R, t, options);
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
    [R, t] = pgo_exp(R, t, x, options);
    fprintf('Iter=%i, cost=%f, gradnorm=%.2e, xnorm=%.2e \n', ...
              iter, cost, gradnorm, norm(x));
end


end