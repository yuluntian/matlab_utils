% function [is_optimal, info] = certify_chordal_rotation_averaging(measurements, Rstar, options)
% Given a chordal rotation averaging problem and a first-order critical
% point Rstar, verify if Rstar is a global minimizer of the optimization problem.
%
% Reference: Eriksson et al. "Rotation Averaging and Strong Duality"
%
% Yulun Tian
function [is_optimal, info] = certify_chordal_rotation_averaging(measurements, Rstar, options)

if nargin < 3
    options = struct;
end
if ~isfield(options, 'symmetric_tol')
    options.symmetric_tol = 1e-4;
end
if ~isfield(options, 'lanczos_tol')
    options.lanczos_tol = 1e-10;
end
if ~isfield(options, 'lanczos_maxiter')
    options.lanczos_maxiter = 1000;
end
if ~isfield(options, 'lanczos_subspace_dim')
    options.lanczos_subspace_dim = 50;
end
if ~isfield(options, 'mineig_tol')
    options.mineig_tol = -1e-8;
end
if ~isfield(options, 'verbose')
    options.verbose = false;
end
edges = measurements.edges;
d = size(measurements.R{1},1);
n = max(max(edges));
m = size(edges, 1);
assert(size(Rstar, 1) == d)
assert(size(Rstar, 2) == d * n)
info = struct;

% Construct data matrix
RTilde = sparse(d*n, d*n);
for k = 1:m
    i = edges(k,1);
    j = edges(k,2);
    assert(i~=j);
    Rij = measurements.R{k};
    check_rotation_matrix(Rij);
    kappa = measurements.kappa{k};
    idxs = (i-1)*d + 1 : i*d;
    jdxs = (j-1)*d + 1 : j*d;
    RTilde(idxs, jdxs) = RTilde(idxs, jdxs) + kappa * Rij;
    RTilde(jdxs, idxs) = RTilde(jdxs, idxs) + kappa * Rij';
end
[RTilde, sym_error] = make_symmetric(RTilde);
assert(sym_error < 1e-6);  % data matrix should always be symmetric

% Recover Langrange multipliers and dual certificate matrix
P = RTilde * Rstar';
S = -RTilde;
multiplier_symmetric_errors = zeros(1,n);
for i = 1:n
    idxs = (i-1)*d + 1 : i*d;
    Pi = P(idxs, :);
    Ri_star = Rstar(:, idxs);
    Lambda_i = Pi * Ri_star;
    [Lambda_i, sym_error] = make_symmetric(Lambda_i);
    multiplier_symmetric_errors(i) = sym_error;
    S(idxs, idxs) = Lambda_i;
end
S = make_symmetric(S);

% Check if dual certificate matrix is positive semidefinite
% TODO: use the spectrum shifting strategy as discussed in SE-Sync?

% Compute the initial guess for eigenvalue computation
% by perturbing a single row of the given first-order critical point
v = Rstar(1,:)';
sigma = 0.2 * norm(v) / (n * d);
vinit = v + sigma*randn(d*n, 1);

[vmin, lambda_min, flag] = eigs(S, 1, 'smallestreal', ...
                                                   'Tolerance', options.lanczos_tol, ...
                                                   'MaxIterations', options.lanczos_maxiter, ...
                                                   'SubspaceDimension', options.lanczos_subspace_dim, ...
                                                   'StartVector', vinit, ...
                                                   'FailureTreatment','keep', ...
                                                   'Display', options.verbose);

info.lambda_min = lambda_min;
info.v_min = vmin;
info.multiplier_symmetric_errors = multiplier_symmetric_errors;
if flag ~=0 
    warning('Eigenvalue computation did not converge to desired precision.')
    is_optimal = false;
elseif max(multiplier_symmetric_errors) > options.symmetric_tol
    [max_err, max_ind] = max(multiplier_symmetric_errors);
    warning('Lagrange multiplier is not sufficiently symmetric. Max error: %.2ef (node %i)', max_err, max_ind)
    is_optimal = false;
elseif lambda_min < options.mineig_tol
    warning('Min eigenvalue is not sufficiently close to zero: %.2e', lambda_min)
    is_optimal = false;
else
    fprintf(['Certified global minimum! ' ...
        'Min eigval: %.2e, ' ...
        'Max multiplier symmetry error: %.2e.' ...
        '\n'], lambda_min, max(multiplier_symmetric_errors));
    is_optimal = true;
end

end