% function [R, info] = rotation_averaging_newton(measurements, R, options)
% A simple implementation of the Newton method for 2D/3D rotation
% averaging. 
% 
% This function supports optimization on the underlying quotient manifold (set
% options.quotient_optimization = true). In this case, the preconditioned conjugate
% gradient (PCG) is used to solve the Newton systems. When quotient
% optimization is disabled, the backslash oprator is used instead, which
% could be numerically unstable.
%
% Yulun Tian
function [R, info] = rotation_averaging_newton(measurements, R, options)
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
if ~isfield(options, 'quotient_optimization')
    options.quotient_optimization = false;
end
if ~isfield(options, 'pcg_tolerance')
    options.pcg_tolerance = 1e-6;
end
if ~isfield(options, 'pcg_maxit')
    options.pcg_maxit = 50;
end

d = size(measurements.R{1}, 1);
n = max(max(measurements.edges));
assert(d == 2 || d == 3, 'Rotation averaging problem must be 2D or 3D.');
if d == 2
    p = 1;
else
    p = 3;
end
info = struct;

if options.quotient_optimization
    fprintf('Performing optimization on quotient manifold.\n');
    options.tangent_space_parametrization = 'global';
    % form orthonormal basis of vertical space
    N = kron(ones(n,1), speye(p));
    for i = 1:p
        N(:,i) = N(:,i) / norm(N(:,i));
    end
    % orthogonal projection operator onto the horizontal space
    Ph_func = @(x) x - N * (N' * x);
    if options.lambda > 0
        warning('Ignoring regularization when optimizing on quotient manifold.');
    end
    % Use regularized laplacian as preconditioner
    if strcmp(options.rotation_distance, 'geodesic')
        weights = [measurements.kappa{:}];
    elseif strcmp(options.rotation_distance, 'chordal')
        weights = 2 * [measurements.kappa{:}];
    else
        error('Unknown rotation distance: %s', options.rotation_distance);
    end
    L = construct_weighted_laplacian(1:n, measurements.edges, weights);
    % Regularization is needed to make sure preconditioner is PD
    precon_lambda = min(1e-8, options.gradnorm_tol);
    Mp = kron(L + precon_lambda * speye(n), speye(p));
end

% Number of Newton steps done.
iter = 0; 
info.cost = [];
info.gradnorm = [];
info.iterations = [];
while true
    cost = evaluate_rotation_averaging_cost(measurements, R, options);
    [grad, Hess] = differentiate_rotation_averaging(measurements, R, options);
    gradnorm = norm(grad);
    info.iterations(end+1) = iter;
    info.cost(end+1) = cost;
    info.gradnorm(end+1) = gradnorm;
    if iter >= options.max_iterations || gradnorm < options.gradnorm_tol
        fprintf('Final result: iter=%i, cost=%f, gradnorm=%.2e. \n', ...
                   iter, cost, gradnorm);
        break;
    end
    
    % Solve Newton system
    % Since the Hessian is always singular (and also dense when using
    % quotient optimization), we will use PCG
    if ~options.quotient_optimization
        % Ignore the quotient structure and perform Newton's method in the total space
        % Without regularization, this could be numerically unstable! 
        if options.lambda == 0
            warning('Solving Newton system using backslash could be unstable!');
        end
        H = Hess + options.lambda * speye(p * n);
        x = - H \ grad;
    else
        % Account for the quotient structure by find the min norm solution to: 
        % (Ph H Ph)x = -grad, 
        % where Ph is the orthogonal projection onto the horizontal space; see [Boumal, eq (9.55)]
        
        % Since the left handside matrix is dense, we will use PCG to avoid
        % explicitly forming the quotient space Hessian
        QH_func = @(x) Ph_func(Hess * Ph_func(x));
        x_pcg = pcg(QH_func, ...
                            -grad, ...
                            options.pcg_tolerance, ...
                            options.pcg_maxit, ...
                            Mp);  

        % Project the solution to the horitontal space
        x = Ph_func(x_pcg);
    end
    % Apply tangent space solution
    R = rotation_averaging_exp(R, x, options);
    fprintf('Iter=%i, cost=%f, gradnorm=%.2e, xnorm=%.2e \n', ...
              iter, cost, gradnorm, norm(x));

    iter = iter + 1;
end

end