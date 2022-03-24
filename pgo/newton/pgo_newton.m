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
assert(d == 3);
p = 6;

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
    for i = 1:n
        idxs = (i-1) * d + 1  : i * d;
        ipxs = (i-1) * p + 1  : i * p;
        Ri = R(:, idxs);
        ti = t(:,i);
        xi = x(ipxs);
        dRi = xi(1:3);
        dti = xi(4:6);
        Ri_new = Ri * exp3(dRi);
        ti_new = ti + dti;
        R(:, idxs) = Ri_new;
        t(:, i) = ti_new;
    end
    fprintf('Iter=%i, cost=%f, gradnorm=%.2e, xnorm=%.2e \n', ...
              iter, cost, gradnorm, norm(x));
end


end