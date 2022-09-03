% Solve the system of linear equation Ax = b using SOR
% Reference: Golub section 11.2
% Input:
% - [A, b]: specify the linear system Ax = b
% - x0: initial estimate
% - options.blocks: blocks of indices
% - options.w: over-relaxation parameter
% - options.maxiter: maximum number of iterations
% - options.gradnorm_tol: gradient norm tolerance
function [x, info] = SOR(A, b, x0, options)
n = size(A,1);
%% Parse options
if ~isfield(options, 'blocks')
    error('Must specify blocks (partitioning of indices) to use SOR!');
else
    blocks = options.blocks;
    k = length(blocks);
end

% Over relaxation parameter
if ~isfield(options, 'w')
    % default to Gauss-Seidel (recommended by Choudhary et al)
    w = 1.0;
else
    w = options.w;
end

% Maximum iterations
if ~isfield(options, 'maxiter')
    maxiter = 1000;
else
    maxiter = options.maxiter;
end

% Stopping condition
if ~isfield(options, 'relchange_tol')
    % Default in Choudhary et al.
    eta = 0.1; 
else
    eta = options.relchange_tol;
end

% flagged initialization
if ~isfield(options, 'flag_init')
    options.flag_init = false;
end
if options.flag_init
    warning('Using flagged initialization will overwrite SOR initial guess.');
    pause;
    x0 = SOR_flagged_initialization(A, b, blocks);
end

% Verbose flag
if ~isfield(options, 'verbose')
    verbose = true;
else
    verbose = options.verbose;
end

% Check that A is symmetric positive definite
if ~isfield(options, 'check_spd')
    check_spd = false;
else
    check_spd = options.check_spd;
end

% Check that splitting is valid (i.e., spectral radius < 1)
if ~isfield(options, 'check_spectral_radius')
    check_spectral_radius = false;
else
    check_spectral_radius = options.check_spectral_radius;
end

% Record iterate
if ~isfield(options, 'save_iterates')
    options.save_iterates = false;
end


%% Sanity checks
% check inputs have consistent shapes
assert(size(A,2) == n);
assert(size(b,1) == n);
assert(size(b,2) == 1);
assert(size(x0,1) == n);
assert(size(x0,1) == n);

% check that A is symmetric positive definite
if check_spd
    [~,p] = chol(A);
    assert(p==0, 'A is not symmetric positive definite!');
end

% check blocks are mutually exclusive
for i = 1:k
    for j = i+1:k
        block_i = blocks{i};
        block_j = blocks{j};
        assert(isempty(intersect(block_i, block_j)), ...
               'input contains overlapping blocks!');
    end
end

% check union of all blocks give full indices
block_all = [];
for i = 1:k
    block_all = union(block_all, blocks{i});
end
block_all = sort(block_all, 'ascend');
assert( all(reshape(block_all, 1, []) == 1:n));

% check each block consists of consecutive indices
for i = 1:k
    assert(all(diff(blocks{i}) == 1),  'Block is not consecutive!');
end

gradnorm = @(x) norm(A*x - b);
% Reference solution
xopt = A \ b;
%% Construct splitting of A

D = sparse(n, n);
L = sparse(n, n);
U = sparse(n, n);
for i = 1:k % block row
   for j = 1:k % block column
       blocki = blocks{i};
       blockj = blocks{j};
       if size(blocki,2) == 0 || size(blockj,2) == 0
           warning('SOR: empty blocks ignored.\n');
           continue;
       end
       
       % Y.T. assume that block consists of consecutive indices
       istart = blocki(1);
       iend = blocki(end);
       jstart = blockj(1);
       jend = blockj(end);
       
       Aij = A(istart:iend, jstart:jend);
       if i < j
           U(istart:iend, jstart:jend) = Aij;
       elseif i > j
           L(istart:iend, jstart:jend) = Aij;
       else
           D(istart:iend, jstart:jend) = Aij;
       end
   end
end
assert(norm(A-U-L-D, 'fro') < 1e-6);
M = (1/w)*D + L;
N = (1/w - 1)*D - U;
% check that M and N is a valid splitting
assert(norm(A - (M-N), 'fro') < 1e-6);

% check spectral radius
if check_spectral_radius
    Gfun = @(x) M \ N * x;
    rho = eigs(Gfun, n, 1, 'largestabs', ...
               'Tolerance', 1e-8, ...
               'MaxIterations', 1000, ...
               'FailureTreatment', 'keep');
    rho = abs(rho);
    fprintf('Spectral radius = %.3f\n', rho);
    assert(rho < 1, 'Splitting is not valid! Spectral radius = %g\n', rho);
end

if verbose
    fprintf('SOR | Finished preprocessing.\n');
end
%% SOR iterations
x = x0;
iterations = [];
distances = [];
residuals = [];
if options.save_iterates
    iterates = [];
end
iter = 0;
while true
    iterations(end+1) = iter;  % count number of SOR steps performed
    distances(end+1) = norm(xopt-x);
    residuals(end+1) = gradnorm(x);
    if options.save_iterates
        iterates = [iterates x];
    end
    
    % stopping condition: relative change
    if iter > 0 && norm(x_prev-x) < eta
        if verbose
            fprintf('SOR | iter = %i | residual = %.5g | rel_change = %.5g | dist_to_opt = %.5g\n', ...
                      iter, residuals(end), norm(x_prev-x), distances(end));
            fprintf('Reached stopping condition.\n');
        end
        break; 
    end
    % stopping condition: maximum iterations
    if iter >= maxiter && verbose
       fprintf('SOR | iter = %i | residual = %.5g | dist_to_opt = %.5g\n', ...
                      iter, residuals(end), distances(end));
       fprintf('SOR | Reached maximum iteration.\n');
       break;
    end
    
    % SOR step
    x_prev = x;
    x = M \ (N*x + b);

    % print progress
    if verbose
        if iter > 0 && mod(iter, 10) == 0
           rel_change = norm(x_prev-x);
           fprintf('SOR | iter = %i | residual = %.5g | rel_change = %.5g | dist_to_opt = %.5g\n', ...
           iter, residuals(end), rel_change, distances(end));
        end
    end
    
    

    iter = iter + 1;
end

info.iterations = iterations;
info.distances = distances;
info.residuals = residuals;
if options.save_iterates
    info.iterates = iterates;
end

end