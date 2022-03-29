% function [R,t] = chordal_pgo_gauss_newton_2d(measurements, R, t, options)
% A simple Gauss-Newton solver for 2D pose graph optimization.
% 
% Yulun Tian
function [R, t] = pgo_gauss_newton(measurements, R, t, options)
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
d = length(measurements.t{1});
n = max(max(measurements.edges));
if d == 2
    p = 3;  % state dimension
    Rexp = @exp2;
elseif d == 3
    p = 6; 
    Rexp = @exp3;
else
    error('Invalid dimension %i', d);
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
    for i = 1:n
        idxs = ((i-1)*d+1) : i*d;
        idxs_ = (i-1) * p + 1: i * p;
        Ri = R(:, idxs);
        ti = t(:,i);
        xi = x(idxs_);
        if d == 2
            dRi = xi(1);
            dti = xi(2:3);
        else
            dRi = xi(1:3);
            dti = xi(4:6);
        end
        ti_new = ti + dti;
        if strcmp(options.tangent_space_parametrization, 'global')
            Ri_new = Rexp(dRi) * Ri;
        elseif strcmp(options.tangent_space_parametrization, 'local')
            Ri_new = Ri * Rexp(dRi);
        else
            error('Unknown tangent space parametrization: %s', options.tangent_space_parametrization);
        end
        R(:, idxs) = Ri_new;
        t(:, i) = ti_new;
    end
    fprintf('Iter=%i, cost=%f, gradnorm=%.2e, xnorm=%.2e \n', ...
              iter, cost, gradnorm, norm(x));
    iter = iter + 1;
end

fprintf('Final result: iter=%i, cost=%f, gradnorm=%.2e. \n', ...
                   iter, cost, gradnorm);

end